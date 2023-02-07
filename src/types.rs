//Various DNA lookup tables and hashing methods are taken from miniprot by Heng Li. Attached below is their license:
//The MIT License

// **** miniprot LICENSE ***
//Copyright (c) 2022-     Dana-Farber Cancer Institute
//
//Permission is hereby granted, free of charge, to any person obtaining
//a copy of this software and associated documentation files (the
//"Software"), to deal in the Software without restriction, including
//without limitation the rights to use, copy, modify, merge, publish,
//distribute, sublicense, and/or sell copies of the Software, and to
//permit persons to whom the Software is furnished to do so, subject to
//the following conditions:
//
//The above copyright notice and this permission notice shall be
//included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//******************************
pub const  DNA_TO_AA: [u8; 64] =
            *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

use std::cmp::Ordering;
use serde::{Deserialize, Serialize};
// bytecheck can be used to validate your data if you want
use smallvec::SmallVec;
use crate::params::*;
use partitions::*;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::str;
pub const BYTE_TO_SEQ: [MarkerBits; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];


pub type GnPosition = u32;
pub type ContigIndex = u32;
//pub type KmerBits = u128;
pub type MarkerBits = u64;
pub type SeedBits = u32;
pub type KmerToSketch = MMHashMap<MarkerBits, SmallVec<[u32; KMER_SK_SMALL_VEC_SIZE]>>;
//pub type KmerToSketch = MMHashMap<MarkerBits, Vec<usize>>;
pub type KmerSeeds = MMHashMap32<SeedBits, SmallVec<[SeedPosition;SMALL_VEC_SIZE]>>;
//pub type KmerSeeds = MMHashMap<SeedBits, SmallVec<[SeedPosition;SMALL_VEC_SIZE]>>;


//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMBuildHasher32 = BuildHasherDefault<MMHasher32>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type MMHashMap32<K, V> = HashMap<K, V, MMBuildHasher32>;
pub type MMHashSet<K> = HashSet<K, MMBuildHasher>;

//Thomas Wang's hash function taken from minimap2
#[inline]
pub fn mm_hashi64(kmer: i64) -> i64 {
    let mut key = kmer as u64;
    key = !(key.wrapping_add(key << 21)); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key as i64
}

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
    key
}

#[inline]
pub fn mm_hash_bytes_32(bytes: &[u8]) -> usize {
    let mut key = (u32::from_ne_bytes(bytes.try_into().unwrap())) as usize;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_ne_bytes(bytes.try_into().unwrap());
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Default, Clone, Serialize, Deserialize, Debug)]
pub struct SeedPosition{
    pub pos: GnPosition,
    pub canonical: bool,
    pub contig_index: ContigIndex,
    pub phase: u8
}

//impl Hash for SeedPosition{
//    fn hash<H: Hasher>(&self, state: &mut H) {
//        self.pos.hash(state);
//    }
//}

#[derive(Eq, PartialEq, Serialize, Deserialize, Clone)]
pub struct Sketch {
    pub file_name: String,
    pub kmer_seeds_k: Option<KmerSeeds>,
    pub contigs: Vec<String>,
    pub total_sequence_length: usize,
    pub contig_lengths: Vec<GnPosition>,

    //repetitive_kmers isn't being used right now, but we may
    //make the repeat masking procedure smarter in the future. 
    //A repetitive k-mer threshold is set in chaining as
    //2500/c. 
    //
    //Also leaving this here so that we users can have 
    //valid .sketch files across versions, although we may
    //break that in the future. 
    pub repetitive_kmers: usize,
    pub marker_seeds: MMHashSet<MarkerBits>,
    pub marker_c: usize,
    pub c: usize,
    pub k: usize,
    pub contig_order: usize,
    pub amino_acid: bool
}

impl Sketch{
    pub fn get_markers_only(sketch: &Sketch) -> Sketch{
        
        Sketch{
            file_name : sketch.file_name.clone(),
            kmer_seeds_k : None,
            contigs: sketch.contigs.clone(),
            total_sequence_length : sketch.total_sequence_length,
            contig_lengths : vec![],
            repetitive_kmers : sketch.repetitive_kmers,
            marker_seeds : sketch.marker_seeds.clone(),
            marker_c : sketch.marker_c,
            c : sketch.c,
            k : sketch.k,
            contig_order: sketch.contig_order,
            amino_acid: sketch.amino_acid
        }
    }
    
    pub fn new(marker_c: usize, c: usize, k: usize, file_name: String, amino_acid: bool) -> Sketch{
        assert!(marker_c >= c);
        let mut new_sketch = Sketch::default();
        new_sketch.c = c;
        new_sketch.k = k;
        new_sketch.marker_c = c;
        new_sketch.file_name = file_name;
        new_sketch.amino_acid = amino_acid;
        new_sketch
    }
}

impl PartialOrd for Sketch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Sketch {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.file_name, self.contig_order).cmp(&(&other.file_name, other.contig_order))
    }
}


impl Hash for Sketch{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.file_name.hash(state);
    }
}
impl Default for Sketch {
    fn default() -> Self {
        Sketch {
            file_name: String::new(),
            kmer_seeds_k: None,
            contigs: vec![],
            total_sequence_length: 0,
            contig_lengths: vec![],
            repetitive_kmers: 0,
            marker_seeds: MMHashSet::default(),
            marker_c: 0,
            c: 0,
            k: 0,
            contig_order:0,
            amino_acid: false,
        }
    }
}

#[derive(Default)]
pub struct MMHasher32 {
    hash: usize,
}

impl Hasher for MMHasher32 {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash_bytes_32(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}





#[derive(Default)]
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



#[derive(Debug, Eq, Hash, Clone)]
pub struct KmerEnc {
    pub kmer: u64,
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
    pub fn decode(byte: u64) -> u8 {
        if byte == 0 {
            b'A'
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

    pub fn print_string(kmer: u64, k: usize) {
        let mut bytes = vec![];
        let mask = 3;
        for i in 0..k {
            let val = kmer >> (2 * i);
            let val = val & mask;
            bytes.push(KmerEnc::decode(val));
        }
        dbg!(str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
    }

    pub fn print_string_aa(kmer: u64, k: usize, sketch_params: &SketchParams) {
        let mut bytes = vec![];
        let mask = 63;
        for i in 0..k {
            let val = kmer >> (6 * i);
            let val = val & mask;
//            bytes.push(KmerEnc::decode(val));
            bytes.push(sketch_params.acgt_to_aa_letters[val as usize]);
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

#[derive(PartialEq, PartialOrd, Debug, Clone, Default)]
pub struct ChainInterval {
    pub score: f64,
    pub num_anchors: usize,
    pub interval_on_query: (GnPosition, GnPosition),
    pub interval_on_ref: (GnPosition, GnPosition),
    pub ref_contig: usize,
    pub query_contig: usize,
    pub chunk_id: usize,
    pub reverse_chain: bool,
    pub overlap: u32
}
impl ChainInterval {
    pub fn query_range_len(&self) -> GnPosition {
        self.interval_on_query.1 - self.interval_on_query.0
    }
    pub fn ref_range_len(&self) -> GnPosition {
        self.interval_on_ref.1 - self.interval_on_ref.0
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
    pub seeds_in_chunk: Vec<Vec<GnPosition>>,
}

#[derive(Default, Clone, PartialEq, Eq, Hash, Debug)]
pub struct Orf{
    pub start: usize,
    pub end: usize,
    pub phase: u8
}

#[derive(Default, Clone, Debug)]
pub struct AniEstResult{
    pub ani: f32,
    pub align_fraction_query: f32,
    pub align_fraction_ref: f32,
    pub ref_file: String,
    pub query_file: String,
    pub query_contig: String,
    pub ref_contig: String,
    pub ci_upper: f32,
    pub ci_lower: f32,
    pub aai: bool,
    pub quant_90_contig_len_q: f32,
    pub quant_90_contig_len_r: f32,
    pub quant_50_contig_len_q: f32,
    pub quant_50_contig_len_r: f32, 
    pub quant_10_contig_len_q: f32,
    pub quant_10_contig_len_r: f32, 
    pub std: f32,
    pub num_contigs_q: u32,
    pub num_contigs_r: u32,
    pub avg_chain_int_len: u32,
    pub total_bases_covered: u32, 
}
