use crate::nthash::NtHashIterator;
use crate::types::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;
use debruijn::Kmer;
use debruijn::Mer;
use debruijn::Vmer;
use fxhash::{hash, FxHashSet};
use std::hash::Hasher;

pub fn fmh_seeds(
    string: &[u8],
    k: usize,
    c: usize,
    contig_index: ContigIndex,
    kmer_seeds: &mut KmerSeeds,
) -> usize{
    let mut num_seeds = 0;
    let mut num_kmers = 0;
    let f_seq_dna_string = DnaString::from_acgt_bytes(string);
    //let r_seq_dna_string = f_seq_dna_string.rc();
    for (i, kmer_f) in f_seq_dna_string.iter_kmers::<VarIntKmer<u64, K30>>().enumerate() {
        num_kmers += 1;
        let (kmer_canonical, canonical) = kmer_f.min_rc_flip();
        let hash = hash(&kmer_canonical);
        if hash < usize::MAX / c {
            let kmer = CanonicalKmer {
                kmer: kmer_canonical,
                canonical: canonical,
            };
            let kmer_positions = kmer_seeds.entry(kmer).or_insert(FxHashSet::default());
            kmer_positions.insert((i.try_into().unwrap(), canonical, contig_index));
            num_seeds += 1;
        }
    }
    dbg!(num_seeds, num_kmers, num_kmers / num_seeds);
    return num_kmers / num_seeds;
}
