//! ntHash is a hash function tuned for genomic data.
//! It performs best when calculating hash values for adjacent k-mers in
//! an input sequence, operating an order of magnitude faster than the best
//! performing alternatives in typical use cases.
//!
//! [Scientific article with more details](https://doi.org/10.1093/bioinformatics/btw397)
//!
//! [Original implementation in C++](https://github.com/bcgsc/ntHash/)
//!
//! This crate is based on ntHash [1.0.4](https://github.com/bcgsc/ntHash/releases/tag/v1.0.4).
//!


pub use crate::error::{Error, Result};

pub(crate) const MAXIMUM_K_SIZE: usize = u32::max_value() as usize;

const H_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'N' as usize] = 0;
    lookup
};

const RC_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'C' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'G' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'T' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'N' as usize] = 0;
    lookup
};

#[inline(always)]
fn h(c: u8) -> u64 {
    let val = H_LOOKUP[c as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", c as char)
    }
    val
}

#[inline(always)]
fn rc(nt: u8) -> u64 {
    let val = RC_LOOKUP[nt as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", nt as char)
    }
    val
}

/// Calculate the hash for a k-mer in the forward strand of a sequence.
///
/// This is a low level function, more useful for debugging than for direct use.
///
/// ```
///    use nthash::ntf64;
///    let fh = ntf64(b"TGCAG", 0, 5);
///    assert_eq!(fh, 0xbafa6728fc6dabf);
/// ```
pub fn ntf64(s: &[u8], i: usize, k: usize) -> u64 {
    let mut out = h(s[i + k - 1]);
    for (idx, v) in s.iter().skip(i).take(k - 1).enumerate() {
        out ^= h(*v).rotate_left((k - idx - 1) as u32);
    }
    out
}

/// Calculate the hash for a k-mer in the reverse strand of a sequence.
///
/// This is a low level function, more useful for debugging than for direct use.
///
/// ```
///    use nthash::ntr64;
///    let rh = ntr64(b"TGCAG", 0, 5);
///    assert_eq!(rh, 0x8cf2d4072cca480e);
/// ```
pub fn ntr64(s: &[u8], i: usize, k: usize) -> u64 {
    let mut out = rc(s[i]);
    for (idx, v) in s.iter().skip(i + 1).take(k - 1).enumerate() {
        out ^= rc(*v).rotate_left(idx as u32 + 1);
    }
    out
}

/// Calculate the canonical hash (minimum hash value between the forward
/// and reverse strands in a sequence).
///
/// This is a low level function, more useful for debugging than for direct use.
///
/// ```
///    use nthash::ntc64;
///    let hash = ntc64(b"TGCAG", 0, 5);
///    assert_eq!(hash, 0xbafa6728fc6dabf);
/// ```
pub fn ntc64(s: &[u8], i: usize, ksize: usize) -> u64 {
    u64::min(ntr64(s, i, ksize), ntf64(s, i, ksize))
}

/// Takes a sequence and ksize and returns the canonical hashes for each k-mer
/// in a Vec. This doesn't benefit from the rolling hash properties of ntHash,
/// serving more for correctness check for the NtHashIterator.
pub fn nthash(seq: &[u8], ksize: usize) -> Vec<u64> {
    seq.windows(ksize).map(|x| ntc64(x, 0, ksize)).collect()
}

/// An efficient iterator for calculating hashes for genomic sequences.
///
/// Since it implements the `Iterator` trait it also
/// exposes many other useful methods. In this example we use `collect` to
/// generate all hashes and put them in a `Vec<u64>`.
/// ```
///     # use nthash::Result;
///     use nthash::NtHashIterator;
///
///     # fn main() -> Result<()> {
///     let seq = b"ACTGC";
///     let iter = NtHashIterator::new(seq, 3)?;
///     let hashes: Vec<u64> = iter.collect();
///     assert_eq!(hashes,
///                vec![0x9b1eda9a185413ce, 0x9f6acfa2235b86fc, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
/// or, in one line:
/// ```
///     # use nthash::Result;
///     use nthash::NtHashIterator;
///
///     # fn main() -> Result<()> {
///     assert_eq!(NtHashIterator::new(b"ACTGC", 3)?.collect::<Vec<u64>>(),
///                vec![0x9b1eda9a185413ce, 0x9f6acfa2235b86fc, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
#[derive(Debug)]
pub struct NtHashIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,
    rh: u64,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> NtHashIterator<'a> {
    /// Creates a new NtHashIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize) -> Result<NtHashIterator<'a>> {
        if k > seq.len() {
            return Err(Error::KSizeOutOfRange {
                ksize: k,
                seq_size: seq.len(),
            });
        }
        if k > MAXIMUM_K_SIZE {
            return Err(Error::KSizeTooBig(k));
        }
        let mut fh = 0;
        for (i, v) in seq[0..k].iter().enumerate() {
            fh ^= h(*v).rotate_left((k - i - 1) as u32);
        }

        let mut rh = 0;
        for (i, v) in seq[0..k].iter().rev().enumerate() {
            rh ^= rc(*v).rotate_left((k - i - 1) as u32);
        }

        Ok(NtHashIterator {
            seq,
            k,
            fh,
            rh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
    }
}

impl<'a> Iterator for NtHashIterator<'a> {
    type Item = (u64,bool);

    fn next(&mut self) -> Option<(u64,bool)> {
        if self.current_idx == self.max_idx {
            return None;
        };

        if self.current_idx != 0 {
            let i = self.current_idx - 1;
            let seqi = self.seq[i];
            let seqk = self.seq[i + self.k];

            self.fh = self.fh.rotate_left(1) ^ h(seqi).rotate_left(self.k as u32) ^ h(seqk);

            self.rh = self.rh.rotate_right(1)
                ^ rc(seqi).rotate_right(1)
                ^ rc(seqk).rotate_left(self.k as u32 - 1);
        }

        self.current_idx += 1;

        let canonical = self.fh > self.rh;
        Some((u64::min(self.rh, self.fh), canonical))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}

impl<'a> ExactSizeIterator for NtHashIterator<'a> {}

/// An efficient iterator for calculating hashes for genomic sequences. This
/// returns the forward hashes, not the canonical hashes.
///
/// Since it implements the `Iterator` trait it also
/// exposes many other useful methods. In this example we use `collect` to
/// generate all hashes and put them in a `Vec<u64>`.
/// ```
///     # use nthash::Result;
///     use nthash::NtHashForwardIterator;
///
///     # fn main() -> Result<()> {
///     let seq = b"ACTGC";
///     let iter = NtHashForwardIterator::new(seq, 3)?;
///     let hashes: Vec<u64> = iter.collect();
///     assert_eq!(hashes, [0xb85d2431d9ba031e, 0xb4d7ab2f9f1306b8, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
/// or, in one line:
/// ```
///     # use nthash::Result;
///     use nthash::NtHashForwardIterator;
///
///     # fn main() -> Result<()> {
///     assert_eq!(NtHashForwardIterator::new(b"ACTGC", 3)?.collect::<Vec<u64>>(),
///                [0xb85d2431d9ba031e, 0xb4d7ab2f9f1306b8, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
#[derive(Debug)]
pub struct NtHashForwardIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> NtHashForwardIterator<'a> {
    /// Creates a new NtHashForwardIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize) -> Result<NtHashForwardIterator<'a>> {
        if k > seq.len() {
            return Err(Error::KSizeOutOfRange {
                ksize: k,
                seq_size: seq.len(),
            });
        }
        if k > MAXIMUM_K_SIZE {
            return Err(Error::KSizeTooBig(k));
        }

        let mut fh = 0;
        for (i, v) in seq[0..k].iter().enumerate() {
            fh ^= h(*v).rotate_left((k - i - 1) as u32);
        }

        Ok(NtHashForwardIterator {
            seq,
            k,
            fh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
    }
}

impl<'a> Iterator for NtHashForwardIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        if self.current_idx == self.max_idx {
            return None;
        };

        if self.current_idx != 0 {
            let i = self.current_idx - 1;
            let seqi = self.seq[i];
            let seqk = self.seq[i + self.k];

            self.fh = self.fh.rotate_left(1) ^ h(seqi).rotate_left(self.k as u32) ^ h(seqk);
        }

        self.current_idx += 1;
        Some(self.fh)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}

impl<'a> ExactSizeIterator for NtHashForwardIterator<'a> {}
