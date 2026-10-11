//! AVX2 forward scan over a sorted seed array.
//!
//! The merge-join in `chain::get_anchors` spends its time answering one
//! question over and over: *given a cursor into the reference's sorted seed
//! array, where is the first seed `>= target`?*
//!
//! For short distances the best answer is a linear scan, not a search. A scan
//! is a sequential stream the hardware prefetcher handles perfectly, whereas
//! an exponential ("galloping") probe breaks that stream and turns each step
//! into an unpredictable access. Vectorising the scan widens the range over
//! which that stays true: one `__m256i` compares 8 seeds, and since the seeds
//! live in their own dense `u32` array (16 per cache line), two vector
//! compares cover a whole cache line.
//!
//! This is only possible because `SeedIndex` keeps seeds separate from
//! positions. Over an interleaved 12-byte record the seeds would be strided
//! and would need a gather, which on this class of hardware costs more than
//! the scalar loads it would replace.

use crate::types::SeedBits;
use std::arch::x86_64::*;

/// Index of the first seed in `seeds[lo..hi]` that is `>= target`, or `hi` if
/// there is none.
///
/// `seeds` must be sorted ascending, though the scan itself only relies on
/// that for the caller's bracket arithmetic.
///
/// # Safety
///
/// Requires AVX2. `lo <= hi <= seeds.len()`.
#[target_feature(enable = "avx2")]
pub unsafe fn scan_forward_avx2(
    seeds: &[SeedBits],
    lo: usize,
    hi: usize,
    target: SeedBits,
) -> usize {
    debug_assert!(lo <= hi && hi <= seeds.len());

    //AVX2 integer comparison is signed, and seeds are u32. Flipping the sign
    //bit maps the unsigned order onto the signed one, so `cmpgt_epi32` then
    //answers the unsigned question.
    let bias = _mm256_set1_epi32(i32::MIN);
    let target_b = _mm256_xor_si256(_mm256_set1_epi32(target as i32), bias);

    let mut i = lo;
    while i + 8 <= hi {
        let v = _mm256_loadu_si256(seeds.as_ptr().add(i) as *const __m256i);
        let v_b = _mm256_xor_si256(v, bias);
        //Lane set <=> seeds[i + lane] < target.
        let lt = _mm256_cmpgt_epi32(target_b, v_b);
        let mask = _mm256_movemask_ps(_mm256_castsi256_ps(lt)) as u32;
        if mask != 0xFF {
            //The first clear bit is the first lane that is >= target.
            return i + mask.trailing_ones() as usize;
        }
        i += 8;
    }

    while i < hi {
        if seeds[i] >= target {
            return i;
        }
        i += 1;
    }
    hi
}
