use debruijn::kmer::*;
use crate::params::*;
use fxhash::{FxHashMap, FxHashSet};
use partitions::*;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::str;
use nohash_hasher::NoHashHasher;




pub const BYTE_TO_SEQ: [KmerBits; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

pub type KmerLength = usize;
pub type KSize = K20;
pub type GnPosition = u32;
pub type ContigIndex = u32;
//pub type KmerBits = u128;
pub type KmerBits = u64;
pub type KmerToSketch = MMHashMap<KmerBits, Vec<Sketch>>;
pub type KmerSeeds = FxHashMap<KmerBits, FxHashSet<SeedPosition>>;


//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type NoHashMap<K, V> = HashMap<K, V, BuildHasherDefault<NoHashHasher<KmerBits>>>;
pub type NoHashSet<K> = HashSet<K, BuildHasherDefault<NoHashHasher<KmerBits>>>;

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_be_bytes(bytes.try_into().unwrap()) as usize;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Default, Clone)]
pub struct SeedPosition{
    pub pos: GnPosition,
    pub canonical: bool,
    pub contig_index: ContigIndex,
    pub phase: u8
}

impl Hash for SeedPosition{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.pos.hash(state);
    }
}

pub struct Sketch {
    pub file_name: String,
    pub kmer_seeds_k: Vec<KmerSeeds>,
    pub contigs: Vec<String>,
    pub total_sequence_length: usize,
    pub repetitive_kmers: Vec<usize>
}
impl Default for Sketch {
    fn default() -> Self {
        return Sketch {
            file_name: String::new(),
            kmer_seeds_k: vec![KmerSeeds::default(); std::mem::size_of::<KmerBits>() * 4],
            contigs: vec![],
            total_sequence_length: 0,
            repetitive_kmers: vec![usize::MAX; std::mem::size_of::<KmerBits>() * 4]
        };
    }
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

#[derive(Debug, Eq, Hash, Clone)]
pub struct KmerEnc {
    pub kmer: KmerBits,
}

//impl Hash for KmerEnc{
//    fn hash<H: Hasher>(&self, state: &mut H) {
//        self.kmer.hash(state);
//    }
//}

impl PartialEq for KmerEnc {
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer
    }
}

impl KmerEnc {
    #[inline]
    pub fn decode(byte: KmerBits) -> u8 {
        if byte == 0 {
            return b'A';
        } else if byte == 1 {
            return b'C';
        } else if byte == 2 {
            return b'G';
        } else if byte == 3 {
            return b'T';
        } else {
            panic!("decoding failed")
        }
    }

    pub fn print_string(kmer: KmerBits, k: usize) {
        let mut bytes = vec![];
        let mask = 3;
        for i in 0..k {
            let val = kmer >> 2 * i;
            let val = val & mask;
            bytes.push(KmerEnc::decode(val));
        }
        dbg!(str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
    }

    pub fn print_string_aa(kmer: KmerBits, k: usize, sketch_params: &SketchParams) {
        let mut bytes = vec![];
        let mask = 63;
        for i in 0..k {
            let val = kmer >> 6 * i;
            let val = val & mask;
//            bytes.push(KmerEnc::decode(val));
            bytes.push(sketch_params.acgt_to_aa_letters[val as usize] as u8);
        }
        dbg!(str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
    }
}

pub struct ChainingResult {
    pub pointer_vec: Vec<usize>,
    pub chain_part: PartitionVec<usize>,
    pub score_vec: Vec<f64>,
    pub num_chunks: usize,
}

pub struct ChainingResultANI {
    pub pointer_vec: Vec<usize>,
}

#[derive(Eq, PartialEq, PartialOrd, Ord, Debug, Default)]
pub struct Anchor {
    pub query_contig: ContigIndex,
    pub query_pos: GnPosition,
    pub ref_contig: ContigIndex,
    pub ref_pos: GnPosition,
    pub ref_phase: u8,
    pub query_phase: u8,
    pub reverse_match: bool,
}

#[derive(PartialEq, PartialOrd, Debug, Clone)]
pub struct ChainInterval {
    pub score: f64,
    pub num_anchors: usize,
    pub interval_on_query: (GnPosition, GnPosition),
    pub interval_on_ref: (GnPosition, GnPosition),
    pub ref_contig: usize,
    pub query_contig: usize,
    pub chunk_id: usize,
}
impl ChainInterval {
    pub fn query_range_len(&self) -> GnPosition {
        return self.interval_on_query.1 - self.interval_on_query.0;
    }
    pub fn ref_range_len(&self) -> GnPosition {
        return self.interval_on_ref.1 - self.interval_on_ref.0;
    }
}

impl Anchor {
    pub fn new(
        rpos: &(GnPosition, ContigIndex),
        qpos: &(GnPosition, ContigIndex),
        ref_phase: u8,
        query_phase: u8,
        reverse: bool,
    ) -> Anchor {
        Anchor {
            ref_pos: rpos.0,
            ref_contig: rpos.1,
            query_pos: qpos.0,
            query_contig: qpos.1,
            ref_phase,
            query_phase,
            reverse_match: reverse,
        }
    }
}

#[derive(Default)]
pub struct AnchorChunks {
    pub chunks: Vec<Vec<Anchor>>,
    pub lengths: Vec<u32>,
    pub seeds_in_chunk: Vec<usize>,
}

#[derive(Default, Clone, PartialEq, Eq, Hash, Debug)]
pub struct Orf{
    pub start: usize,
    pub end: usize,
    pub phase: u8
}

#[derive(Default, Clone, Debug)]
pub struct AniEstResult{
    pub ani: f64,
    pub align_fraction: f64,
    pub ref_file: String,
    pub query_file: String,
    pub query_contig: String,
}