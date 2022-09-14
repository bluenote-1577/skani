use debruijn::kmer::*;
use debruijn::Kmer;
use fxhash::{FxHashSet, FxHashMap};
use std::collections::HashMap;
use std::hash::{Hash, Hasher,BuildHasherDefault};

pub type KSize = K30;
pub type GnPosition = u32;
pub type ContigIndex = u32;
pub type KmerToSketch = FxHashMap<CanonicalKmer, Vec<Sketch>>;
pub type KmerSeeds = FxHashMap<CanonicalKmer, FxHashSet<(GnPosition, bool, ContigIndex)>>;

//Implement minimap2 hashing, will test later. 
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V>  = HashMap<K, V, MMBuildHasher>;


pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_be_bytes(bytes.try_into().unwrap());
    key = !key + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ key >> 28;
    key = key + (key << 31);
    return key;
}

#[derive(Default)]
pub struct Sketch{
    pub file_name: String,
    pub kmer_seeds: KmerSeeds,
    pub contigs: Vec<String>,
    pub total_sequence_length: usize,
    pub c_adj: usize
}

pub struct MMHasher {
    hash: usize,
}

impl Hasher for MMHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl Default for MMHasher {
    #[inline]
    fn default() -> MMHasher {
        MMHasher { hash: 0 }
    }
}


#[derive(Debug)]
pub struct CanonicalKmer {
    pub kmer: VarIntKmer<u64, KSize>,
    pub canonical: bool,
}


impl Hash for CanonicalKmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kmer.hash(state);
    }
}

impl PartialEq for CanonicalKmer {
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer
    }
}
impl Eq for CanonicalKmer {}

#[derive(Eq, PartialEq, PartialOrd, Ord)]
pub struct Anchor{
    pub ref_contig: ContigIndex,
    pub query_contig: ContigIndex,
    pub ref_pos: GnPosition,
    pub query_pos: GnPosition,
    pub reverse_match: bool
}

impl Anchor{
    pub fn new(rpos: &(GnPosition,ContigIndex), qpos: &(GnPosition,ContigIndex), reverse: bool) -> Anchor{
        Anchor{ ref_pos: rpos.0, ref_contig: rpos.1, query_pos: qpos.0, query_contig: qpos.1, reverse_match: reverse}
    }
}
