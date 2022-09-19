use crate::types::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;
use debruijn::Kmer;
use debruijn::Mer;
use debruijn::Vmer;
use fxhash::{hash, FxHashSet};
use std::hash::Hasher;

fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn fmh_seeds(
    string: &[u8],
    k: usize,
    c: usize,
    contig_index: ContigIndex,
    kmer_seeds: &mut KmerSeeds,
) -> f64 {
    if string.len() < 2*k{
        return 1.
    }
    let mut num_seeds = 0;
    let mut rolling_kmer_f: KmerBits = 0;
    let mut rolling_kmer_r: KmerBits = 0;
    let km = k-1;
    let mask = KmerBits::MAX >> (std::mem::size_of::<KmerBits>() * 8 - 2 * k);
    let rc_mask = !(0 | (3 << 2 * (k-1)));
    let len = string.len();
    let threshold = u64::MAX / c as u64;
    for i in 0..len {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
//        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= mask;
//        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= rc_mask;
        rolling_kmer_r |= nuc_r << 2 * (km);
        rolling_kmer_r &= mask;

        if i < k-1{
            continue
        }
//        KmerEnc::print_string(rolling_kmer_f, k);
//        KmerEnc::print_string(rolling_kmer_r, k);
        let canonical = rolling_kmer_f < rolling_kmer_r;
        let canonical_kmer = if canonical{
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };
//        let hash = mm_hash128(canonical_kmer) as KmerBits;
        let hash = mm_hash64(canonical_kmer as u64);
//        let hash = hash(&canonical_kmer) as u64;
//        KmerEnc::print_string(canonical_kmer, k);
        if hash < threshold{
            let kmer = KmerEnc {
                kmer: canonical_kmer,
            };
            num_seeds += 1;
            let kmer_positions = kmer_seeds.entry(kmer).or_insert(FxHashSet::default());
            kmer_positions.insert((i.try_into().unwrap(), canonical, contig_index));
        }
    }
    let c_adj = (string.len() - k + 1) as f64 / num_seeds as f64;
    return c_adj
}

pub fn open_sync_seeds(
    string: &[u8],
    k: usize,
    c: usize,
    contig_index: ContigIndex,
    kmer_seeds: &mut KmerSeeds,
) -> f64 {
    let s = 8;
    let mut num_seeds = 0;
    let w = k - s + 1;
    let t = w / 2 + 1;
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];
    let f_seq_dna_string = DnaString::from_acgt_bytes(string);
    let threshold = usize::MAX / c ;
    //let r_seq_dna_string = f_seq_dna_string.rc();
    //
    for i in 0..f_seq_dna_string.len() - s + 1 {
        let smer: Kmer8 = f_seq_dna_string.slice(i, i + s).get_kmer(0);
        let (hash_smer, _b) = smer.min_rc_flip();
        window_hashes[running_pos] = hash(&hash_smer);
        if i < w - 1 {
            continue;
        }

        if min_running_pos == usize::MAX {
            min_running_pos = position_min(&window_hashes).unwrap();
        } else {
            if min_running_pos == running_pos {
                min_running_pos = position_min(&window_hashes).unwrap();
            } else {
                if window_hashes[running_pos] < window_hashes[min_running_pos] {
                    min_running_pos = running_pos;
                }
            }
        }

        if running_pos > min_running_pos {
            if running_pos - min_running_pos == t - 1 {
                let kmer_f: VarIntKmer<u64, KSize> =
                    f_seq_dna_string.slice(i - w + 1, i - w + 1 + k).get_kmer(0);
                let (kmer_canonical, canonical) = kmer_f.min_rc_flip();
                let downsample;
                if hash(&kmer_canonical) < threshold{
                    downsample = true;
                } else {
                    downsample = false;
                }
                if downsample {
                    let kmer = KmerEnc{
                        kmer: 0,
                    };
                    let kmer_positions = kmer_seeds.entry(kmer).or_insert(FxHashSet::default());
                    kmer_positions.insert((
                        (i + 1 - w).try_into().unwrap(),
                        canonical,
                        contig_index,
                    ));

                    num_seeds += 1;
                }
            }
        } else {
            if w - (min_running_pos - running_pos) == t - 1 {
                let kmer_f: VarIntKmer<u64, KSize> =
                    f_seq_dna_string.slice(i + 1 - w, i + 1 + k - w).get_kmer(0);
                let (kmer_canonical, canonical) = kmer_f.min_rc_flip();
                let downsample;
                if hash(&kmer_canonical) < threshold{
                    downsample = true;
                } else {
                    downsample = false;
                }
                if downsample {
                    let kmer = KmerEnc{
                        kmer: 0,
                    };
                    let kmer_positions = kmer_seeds.entry(kmer).or_insert(FxHashSet::default());
                    kmer_positions.insert((
                        (i + 1 - w).try_into().unwrap(),
                        canonical,
                        contig_index,
                    ));
                    num_seeds += 1;
                }
            }
        }

        running_pos += 1;
        running_pos %= w;
    }

    dbg!(f_seq_dna_string.len() / num_seeds);
    return f_seq_dna_string.len() as f64 / num_seeds as f64;
}
