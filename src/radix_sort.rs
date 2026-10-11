//! LSD radix sorts for the two hot sorts in `chain::get_anchors`.
//!
//! Anchor sorting is the single most expensive thing skani does on a
//! high-identity comparison: at 99.8% ANI it is ~65% of `get_anchors`, which
//! is itself ~64% of runtime. A comparison sort over 27-45k 20-byte structs
//! with a five-field derived `Ord` costs ~400k branchy comparisons; a radix
//! sort over the packed key costs three or four linear passes.
//!
//! Both sorts here reproduce `sort_unstable()`'s result exactly, not merely an
//! equivalent ordering - see `sort_anchors` for why the key alone is not
//! enough and what fixes it.
//!
//! Scratch buffers are thread-local and persist across comparisons, so the
//! steady state allocates nothing.

use crate::params::*;
use crate::types::*;
use std::cell::RefCell;

thread_local! {
    static ANCHOR_SCRATCH: RefCell<Vec<Anchor>> = RefCell::new(Vec::new());
    static U32_SCRATCH: RefCell<Vec<GnPosition>> = RefCell::new(Vec::new());
}

const RADIX_BITS: usize = 8;
const RADIX_BUCKETS: usize = 1 << RADIX_BITS;
const RADIX_MASK: u64 = (RADIX_BUCKETS - 1) as u64;

/// Bits needed to represent `value`, minimum 1.
#[inline]
fn bit_width(value: u32) -> u32 {
    u32::max(32 - value.leading_zeros(), 1)
}

/// Sort anchors into the order `sort_unstable()` would produce.
///
/// The radix key packs only `(query_contig, query_pos)`, the first two fields
/// of `Anchor`'s derived `Ord`. `query_pos` is placed in the low bits using
/// exactly as many as the largest value needs, so the packed key is
/// order-preserving *and* usually fits in three or four bytes - a contig index
/// shifted by a fixed 32 would have forced five passes for no reason.
///
/// Anchors sharing a key are then ordered by `ref_contig`, `ref_pos` and
/// `reverse_match`, which the key does not capture, so each run of equal keys
/// is re-sorted with the full comparison afterwards. Those runs only arise
/// when one query position matches several reference positions, so they are
/// almost always length one and the pass is close to free. Because the full
/// `Ord` is a total order and any two genuinely equal anchors are
/// interchangeable, the result is uniquely determined - this is not merely
/// "a" valid sort order, it is the same one.
pub fn sort_anchors(anchors: &mut Vec<Anchor>) {
    let n = anchors.len();
    if n < RADIX_SORT_MIN_ANCHORS {
        anchors.sort_unstable();
        return;
    }

    let mut max_pos = 0u32;
    let mut max_contig = 0u32;
    for anchor in anchors.iter() {
        max_pos = u32::max(max_pos, anchor.query_pos);
        max_contig = u32::max(max_contig, anchor.query_contig);
    }
    let pos_bits = bit_width(max_pos);
    let key_bits = pos_bits + bit_width(max_contig);
    let passes = ((key_bits as usize) + RADIX_BITS - 1) / RADIX_BITS;

    let key = |anchor: &Anchor| -> u64 {
        ((anchor.query_contig as u64) << pos_bits) | (anchor.query_pos as u64)
    };

    ANCHOR_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        if scratch.len() < n {
            scratch.resize(n, Anchor::default());
        }

        for pass in 0..passes {
            let shift = pass * RADIX_BITS;

            let mut counts = [0u32; RADIX_BUCKETS];
            for anchor in anchors.iter() {
                counts[((key(anchor) >> shift) & RADIX_MASK) as usize] += 1;
            }
            //A pass whose every key lands in one bucket cannot reorder
            //anything; skipping it also skips a full copy of the data.
            if counts.iter().any(|&c| c as usize == n) {
                continue;
            }

            let mut running = 0u32;
            for count in counts.iter_mut() {
                let this = *count;
                *count = running;
                running += this;
            }

            for anchor in anchors.iter() {
                let bucket = ((key(anchor) >> shift) & RADIX_MASK) as usize;
                scratch[counts[bucket] as usize] = *anchor;
                counts[bucket] += 1;
            }
            //Swap the sorted image back into `anchors`; `scratch` keeps its
            //allocation for the next pass and the next comparison.
            anchors[..n].copy_from_slice(&scratch[..n]);
        }
    });

    //Order the fields the key does not cover.
    let mut i = 0;
    while i < n {
        let k = key(&anchors[i]);
        let mut j = i + 1;
        while j < n && key(&anchors[j]) == k {
            j += 1;
        }
        if j - i > 1 {
            anchors[i..j].sort_unstable();
        }
        i = j;
    }
}

/// Sort genome positions ascending, matching `sort_unstable()`.
///
/// Used for the per-contig query position lists, which are 12% of
/// `get_anchors` at high identity and 26% at low identity.
pub fn sort_positions(positions: &mut Vec<GnPosition>) {
    let n = positions.len();
    if n < RADIX_SORT_MIN_POSITIONS {
        positions.sort_unstable();
        return;
    }

    let mut max_value = 0u32;
    for value in positions.iter() {
        max_value = u32::max(max_value, *value);
    }
    let passes = ((bit_width(max_value) as usize) + RADIX_BITS - 1) / RADIX_BITS;

    U32_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        if scratch.len() < n {
            scratch.resize(n, 0);
        }

        for pass in 0..passes {
            let shift = pass * RADIX_BITS;

            let mut counts = [0u32; RADIX_BUCKETS];
            for value in positions.iter() {
                counts[((*value as u64 >> shift) & RADIX_MASK) as usize] += 1;
            }
            if counts.iter().any(|&c| c as usize == n) {
                continue;
            }

            let mut running = 0u32;
            for count in counts.iter_mut() {
                let this = *count;
                *count = running;
                running += this;
            }

            for value in positions.iter() {
                let bucket = ((*value as u64 >> shift) & RADIX_MASK) as usize;
                scratch[counts[bucket] as usize] = *value;
                counts[bucket] += 1;
            }
            positions[..n].copy_from_slice(&scratch[..n]);
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_anchor(contigs: u32, span: u32) -> Anchor {
        Anchor {
            query_contig: fastrand::u32(0..contigs),
            query_pos: fastrand::u32(0..span),
            ref_contig: fastrand::u32(0..contigs),
            ref_pos: fastrand::u32(0..span),
            reverse_match: fastrand::bool(),
        }
    }

    /// The radix sort must agree with `sort_unstable` element for element, not
    /// just as a multiset: `pointer_vec` backtracking and the chunking both
    /// depend on the exact order. Narrow position spans are included on
    /// purpose, since they force the equal-key runs that the key alone cannot
    /// order.
    #[test]
    fn radix_anchor_sort_matches_sort_unstable() {
        fastrand::seed(4242);
        for &n in &[0usize, 1, 2, 3, 255, 256, 257, 1000, 27000] {
            for &(contigs, span) in &[(1u32, 4u32), (1, 1_000_000), (3, 50), (31, 150_000), (1000, 10)]
            {
                let mut a: Vec<Anchor> = (0..n).map(|_| random_anchor(contigs, span)).collect();
                let mut want = a.clone();
                want.sort_unstable();
                sort_anchors(&mut a);
                assert_eq!(a, want, "n={} contigs={} span={}", n, contigs, span);
            }
        }
    }

    #[test]
    fn radix_position_sort_matches_sort_unstable() {
        fastrand::seed(99);
        for &n in &[0usize, 1, 2, 255, 256, 257, 5000, 40000] {
            for &span in &[1u32, 2, 1000, u32::MAX] {
                let mut p: Vec<GnPosition> = (0..n).map(|_| fastrand::u32(0..span.max(1))).collect();
                let mut want = p.clone();
                want.sort_unstable();
                sort_positions(&mut p);
                assert_eq!(p, want, "n={} span={}", n, span);
            }
        }
    }
}
