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
//LEGACY ON-DISK TYPES. The in-memory seed index is `SeedIndex` (a flat, sorted
//vector); these two types only exist so that `.sketch` files, `sketches.db` and
//`markers.bin` keep exactly the same bincode layout as skani <= 0.3.2. They are
//used by `SketchSerde` and converted to/from `SeedIndex` on load/store.
pub type KmerSeeds = MMHashMap32<SeedBits, u64>; // Tagged index: bit 0 = single/multiple, bits 1-63 = data
pub type MultiPositionStorage = Vec<SmallVec<[SeedPosition; 3]>>; // For cases with multiple SeedPositions


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

#[derive(Eq, PartialEq, Ord, PartialOrd, Default, Clone, Copy, Serialize, Deserialize, Debug)]
pub struct SeedPosition{
    pub pos: GnPosition,
    pub contig_index_canonical: ContigIndexCanonical, // packed: contig_index + canonical bit
}

// Type alias for packed contig index with canonical bit in LSB
pub type ContigIndexCanonical = u32;

impl SeedPosition {
    /// Create a new SeedPosition with 30-bit contig index
    pub fn new(pos: GnPosition, contig_index: ContigIndex, canonical: bool) -> Self {
        // Limit contig_index to 30 bits (0x3FFFFFFF)
        assert!(contig_index < (1u32 << 30), "contig_index exceeds 30-bit limit");
        let contig_index_canonical = (contig_index << 1) | (canonical as u32);
        Self {
            pos,
            contig_index_canonical,
        }
    }
    
    /// Extract the canonical flag from the packed field
    #[inline]
    pub fn canonical(&self) -> bool {
        (self.contig_index_canonical & 1) != 0
    }
    
    /// Extract the contig index from the packed field (30 bits)
    #[inline]
    pub fn contig_index(&self) -> ContigIndex {
        self.contig_index_canonical >> 1
    }
    
    /// Set the canonical flag
    #[inline]
    pub fn set_canonical(&mut self, canonical: bool) {
        if canonical {
            self.contig_index_canonical |= 1;
        } else {
            self.contig_index_canonical &= !1;
        }
    }
    
    /// Set the contig index (30 bits max)
    #[inline]
    pub fn set_contig_index(&mut self, contig_index: ContigIndex) {
        debug_assert!(contig_index < (1u32 << 30), "contig_index exceeds 30-bit limit");
        let canonical_bit = self.contig_index_canonical & 1;
        self.contig_index_canonical = (contig_index << 1) | canonical_bit;
    }
    
    /// Pack SeedPosition into 63 bits for tagged index
    #[inline]
    pub fn pack_to_u64(&self) -> u64 {
        // pos: 32 bits, contig_index_canonical: 31 bits (30 + 1 canonical)
        ((self.pos as u64) << 31) | (self.contig_index_canonical as u64)
    }
    
    /// Unpack SeedPosition from 63 bits
    #[inline]
    pub fn unpack_from_u64(packed: u64) -> Self {
        let pos = (packed >> 31) as GnPosition;
        let contig_index_canonical = (packed & 0x7FFFFFFF) as ContigIndexCanonical;
        Self {
            pos,
            contig_index_canonical,
        }
    }
}

/// One (seed, position) record of the flat seed index.
///
/// 12 bytes, and the derived `Ord` orders by `(seed, pos, contig_index_canonical)`,
/// which is exactly the order the index is kept in.
#[derive(Eq, PartialEq, Ord, PartialOrd, Default, Clone, Copy, Debug)]
pub struct SeedEntry {
    pub seed: SeedBits,
    pub position: SeedPosition,
}

/// Flat, sorted, per-genome seed index.
///
/// Replaces the old `HashMap<u32, tagged-index>` + `Vec<SmallVec<[SeedPosition; 3]>>`
/// pair. All `(seed, position)` records live in one contiguous allocation sorted by
/// `(seed, pos, contig)`, so equal seeds form a contiguous run. Anchor generation
/// between two genomes is then a sequential merge-join over two such vectors instead
/// of a hash probe per seed, and there is no per-lookup allocation at all.
///
/// Build order: `push()` while sketching each contig, then `finalize()` once per
/// genome. `finalize()` is what establishes the sorted invariant; everything that
/// reads the index requires a finalized one.
#[derive(Eq, PartialEq, Clone, Default, Debug)]
pub struct SeedIndex {
    //Struct-of-arrays, not an array of `SeedEntry`. The merge-join searches
    //only the seed keys, so keeping them in their own dense `u32` array puts
    //16 of them in a cache line instead of the 5.33 a 12-byte interleaved
    //record manages. That is a 3x reduction in both the footprint being
    //searched and the lines touched per probe, which is what keeps the search
    //in cache for large references, and it is also what makes a vectorised
    //scan possible at all - gathering a strided field would cost more than
    //the scalar loads it replaces.
    //
    //`seeds[i]` and `positions[i]` describe the same record; both are sorted
    //by `(seed, pos, contig)`.
    seeds: Vec<SeedBits>,
    positions: Vec<SeedPosition>,
    //Scratch used only while sketching, emptied by `finalize`.
    build: Vec<SeedEntry>,
    num_distinct_seeds: u32,
    finalized: bool,
}

impl SeedIndex {
    /// Append a record. Cheap; ordering is established later by `finalize`.
    #[inline]
    pub fn push(&mut self, seed: SeedBits, position: SeedPosition) {
        debug_assert!(!self.finalized, "pushing into an already finalized SeedIndex");
        self.build.push(SeedEntry { seed, position });
    }

    /// Sort the records, split them into the two parallel arrays and cache the
    /// distinct-seed count. Idempotent.
    ///
    /// The sort runs over the interleaved `build` buffer so that the ordering
    /// is by the whole record, then the buffer is split and released. That
    /// costs a transient 2x on the index during sketching only; the resident
    /// form is the same 12 bytes per record as the interleaved one.
    pub fn finalize(&mut self) {
        if self.finalized {
            return;
        }
        //Sorting on the whole entry (not just the seed) keeps the index
        //independent of the sort algorithm and of the order in which contigs
        //and SIMD lanes happened to emit their seeds.
        self.build.sort_unstable();
        self.seeds.clear();
        self.positions.clear();
        self.seeds.reserve_exact(self.build.len());
        self.positions.reserve_exact(self.build.len());
        let mut distinct = 0u32;
        let mut prev: Option<SeedBits> = None;
        for e in self.build.iter() {
            if prev != Some(e.seed) {
                distinct += 1;
                prev = Some(e.seed);
            }
            self.seeds.push(e.seed);
            self.positions.push(e.position);
        }
        self.build = Vec::new();
        self.num_distinct_seeds = distinct;
        self.finalized = true;
    }

    /// The sorted seed keys. Only valid after `finalize`.
    #[inline]
    pub fn seeds(&self) -> &[SeedBits] {
        debug_assert!(self.finalized, "SeedIndex read before finalize()");
        &self.seeds
    }

    /// The positions, parallel to `seeds`. Only valid after `finalize`.
    #[inline]
    pub fn positions(&self) -> &[SeedPosition] {
        debug_assert!(self.finalized, "SeedIndex read before finalize()");
        &self.positions
    }

    /// Number of records.
    #[inline]
    pub fn len(&self) -> usize {
        self.seeds.len()
    }

    /// The records grouped by seed, in seed order.
    pub fn runs(&self) -> SeedRuns<'_> {
        debug_assert!(self.finalized, "SeedIndex read before finalize()");
        SeedRuns {
            seeds: &self.seeds,
            positions: &self.positions,
            next: 0,
        }
    }

    /// Number of distinct seeds (the old `kmer_seeds_k.len()`).
    #[inline]
    pub fn num_distinct_seeds(&self) -> usize {
        self.num_distinct_seeds as usize
    }

    /// Build from the legacy tagged-hash-map representation found in `.sketch` files.
    fn from_legacy(kmer_seeds: KmerSeeds, multi_position_storage: &MultiPositionStorage) -> Self {
        let mut index = SeedIndex::default();
        index.build.reserve(kmer_seeds.len());
        for (seed, tagged_index) in kmer_seeds.into_iter() {
            if TaggedIndex::is_single(tagged_index) {
                index.push(seed, TaggedIndex::get_single(tagged_index));
            } else {
                let storage_index = TaggedIndex::get_storage_index(tagged_index);
                for position in multi_position_storage[storage_index].iter() {
                    index.push(seed, *position);
                }
            }
        }
        index.finalize();
        index
    }

    /// Number of distinct seeds that occur more than once, i.e. the length of the
    /// legacy `multi_position_storage` vector this index would render into.
    fn num_multi_seeds(&self) -> usize {
        self.runs().filter(|run| run.positions.len() > 1).count()
    }
}

/// One distinct seed and all of its positions.
pub struct SeedRun<'a> {
    pub seed: SeedBits,
    pub positions: &'a [SeedPosition],
}

/// Iterator over a `SeedIndex`'s records, grouped by seed.
pub struct SeedRuns<'a> {
    seeds: &'a [SeedBits],
    positions: &'a [SeedPosition],
    next: usize,
}

impl<'a> Iterator for SeedRuns<'a> {
    type Item = SeedRun<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.seeds.len() {
            return None;
        }
        let start = self.next;
        let seed = self.seeds[start];
        let mut end = start + 1;
        while end < self.seeds.len() && self.seeds[end] == seed {
            end += 1;
        }
        self.next = end;
        Some(SeedRun {
            seed,
            positions: &self.positions[start..end],
        })
    }
}

/// Truncate contig name to first part (before first whitespace) if short_header is enabled
pub fn truncate_contig_name(name: &str, short_header: bool) -> String {
    if short_header {
        name.split_whitespace().next().unwrap_or(name).to_string()
    } else {
        name.to_string()
    }
}


/// Tagged index utilities for KmerSeeds
pub struct TaggedIndex;

impl TaggedIndex {
    const SINGLE_BIT: u64 = 1;
    
    /// Create tagged index for single SeedPosition
    #[inline]
    pub fn single(seed_position: &SeedPosition) -> u64 {
        Self::SINGLE_BIT | (seed_position.pack_to_u64() << 1)
    }
    
    /// Create tagged index for multiple positions (index into MultiPositionStorage)
    #[inline]
    pub fn multiple(storage_index: usize) -> u64 {
        (storage_index as u64) << 1 // bit 0 = 0 for multiple
    }
    
    /// Check if tagged index represents single SeedPosition
    #[inline]
    pub fn is_single(tagged_index: u64) -> bool {
        (tagged_index & Self::SINGLE_BIT) != 0
    }
    
    /// Extract single SeedPosition from tagged index
    #[inline]
    pub fn get_single(tagged_index: u64) -> SeedPosition {
        debug_assert!(Self::is_single(tagged_index));
        let packed = tagged_index >> 1;
        SeedPosition::unpack_from_u64(packed)
    }
    
    /// Extract storage index for multiple positions
    #[inline]
    pub fn get_storage_index(tagged_index: u64) -> usize {
        debug_assert!(!Self::is_single(tagged_index));
        (tagged_index >> 1) as usize
    }
}

//impl Hash for SeedPosition{
//    fn hash<H: Hasher>(&self, state: &mut H) {
//        self.pos.hash(state);
//    }
//}

#[derive(Eq, PartialEq, Deserialize, Clone)]
#[serde(from = "SketchSerde")]
pub struct Sketch {
    pub file_name: String,
    pub seed_index: Option<SeedIndex>,
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
    pub individual_contig: bool,
    pub amino_acid: bool
}

impl Sketch{
    /// Append a (seed, position) record to this sketch's flat seed index.
    ///
    /// Records are appended unsorted; `finalize_seed_index` must be called once
    /// the genome (or individual contig, in `-i` mode) has been fully sketched.
    #[inline]
    pub fn add_seed_position(&mut self, seed: SeedBits, position: SeedPosition) {
        if let Some(seed_index) = &mut self.seed_index {
            seed_index.push(seed, position);
        }
    }

    /// Sort and freeze the seed index. Must be called after sketching and before
    /// the sketch is used for mapping. Idempotent and a no-op for marker-only
    /// sketches, which carry no seed index.
    pub fn finalize_seed_index(&mut self) {
        if let Some(seed_index) = &mut self.seed_index {
            seed_index.finalize();
        }
    }

    /// Number of distinct seeds in this sketch (0 if it has no seed index).
    pub fn num_distinct_seeds(&self) -> usize {
        self.seed_index
            .as_ref()
            .map_or(0, |index| index.num_distinct_seeds())
    }

    pub fn get_markers_only(sketch: &Sketch) -> Sketch{

        Sketch{
            file_name : sketch.file_name.clone(),
            seed_index : None,
            contigs: sketch.contigs.clone(),
            total_sequence_length : sketch.total_sequence_length,
            contig_lengths : vec![],
            repetitive_kmers : sketch.repetitive_kmers,
            marker_seeds : sketch.marker_seeds.clone(),
            marker_c : sketch.marker_c,
            c : sketch.c,
            k : sketch.k,
            contig_order: sketch.contig_order,
            individual_contig: sketch.individual_contig,
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

/// On-disk image of a `Sketch`, used for reading.
///
/// This mirrors the pre-flat-index `Sketch` layout field for field, in the same
/// order, with the same types, so deserializing it accepts exactly what skani
/// <= 0.3.2 wrote. `Sketch` deserializes through it via `#[serde(from)]` and the
/// flat `SeedIndex` is built from the tagged hash map as the last step of the
/// load. Writing goes the other way through `impl Serialize for Sketch`, which
/// emits this same layout without materializing it.
#[derive(Deserialize)]
struct SketchSerde {
    file_name: String,
    kmer_seeds_k: Option<KmerSeeds>,
    multi_position_storage: MultiPositionStorage,
    contigs: Vec<String>,
    total_sequence_length: usize,
    contig_lengths: Vec<GnPosition>,
    repetitive_kmers: usize,
    marker_seeds: MMHashSet<MarkerBits>,
    marker_c: usize,
    c: usize,
    k: usize,
    contig_order: usize,
    individual_contig: bool,
    amino_acid: bool,
}

impl From<SketchSerde> for Sketch {
    fn from(disk: SketchSerde) -> Sketch {
        let seed_index = disk
            .kmer_seeds_k
            .map(|kmer_seeds| SeedIndex::from_legacy(kmer_seeds, &disk.multi_position_storage));
        Sketch {
            file_name: disk.file_name,
            seed_index,
            contigs: disk.contigs,
            total_sequence_length: disk.total_sequence_length,
            contig_lengths: disk.contig_lengths,
            repetitive_kmers: disk.repetitive_kmers,
            marker_seeds: disk.marker_seeds,
            marker_c: disk.marker_c,
            c: disk.c,
            k: disk.k,
            contig_order: disk.contig_order,
            individual_contig: disk.individual_contig,
            amino_acid: disk.amino_acid,
        }
    }
}

/// Writes the legacy `kmer_seeds_k` field: `Option` of a map from seed to tagged
/// index. Streamed straight out of the sorted index, so nothing is materialized.
/// Storage indices are handed out in seed order, which is the same order
/// `LegacyMultiStorage` emits the position vectors in.
struct LegacyKmerSeeds<'a>(Option<&'a SeedIndex>);

impl Serialize for LegacyKmerSeeds<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self.0 {
            None => serializer.serialize_none(),
            Some(index) => serializer.serialize_some(&LegacyKmerSeedsMap(index)),
        }
    }
}

struct LegacyKmerSeedsMap<'a>(&'a SeedIndex);

impl Serialize for LegacyKmerSeedsMap<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeMap;
        let mut map = serializer.serialize_map(Some(self.0.num_distinct_seeds()))?;
        let mut storage_index = 0usize;
        for run in self.0.runs() {
            let tagged = if run.positions.len() == 1 {
                TaggedIndex::single(&run.positions[0])
            } else {
                let tagged = TaggedIndex::multiple(storage_index);
                storage_index += 1;
                tagged
            };
            map.serialize_entry(&run.seed, &tagged)?;
        }
        map.end()
    }
}

/// Writes the legacy `multi_position_storage` field: the position vectors of the
/// seeds that occur more than once, in the order `LegacyKmerSeeds` numbered them.
struct LegacyMultiStorage<'a>(Option<&'a SeedIndex>);

impl Serialize for LegacyMultiStorage<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let index = match self.0 {
            None => return serializer.serialize_seq(Some(0))?.end(),
            Some(index) => index,
        };
        let mut seq = serializer.serialize_seq(Some(index.num_multi_seeds()))?;
        for run in index.runs().filter(|run| run.positions.len() > 1) {
            seq.serialize_element(&LegacyPositions(run.positions))?;
        }
        seq.end()
    }
}

/// One position vector. `SmallVec` serializes as a plain sequence, so this is
/// byte-compatible with the `SmallVec<[SeedPosition; 3]>` it replaces.
struct LegacyPositions<'a>(&'a [SeedPosition]);

impl Serialize for LegacyPositions<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for position in self.0.iter() {
            seq.serialize_element(position)?;
        }
        seq.end()
    }
}

/// Writes the pre-flat-index `Sketch` layout: same fields, same order, same
/// types as `SketchSerde`, so the byte stream is exactly what skani <= 0.3.2
/// produced and what it can still read back. The seed index is streamed out of
/// its sorted form rather than being rebuilt as a hash map first, which also
/// makes the written bytes deterministic (the old `HashMap` iteration order was
/// not).
impl Serialize for Sketch {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let index = self.seed_index.as_ref();
        let mut sketch = serializer.serialize_struct("Sketch", 14)?;
        sketch.serialize_field("file_name", &self.file_name)?;
        sketch.serialize_field("kmer_seeds_k", &LegacyKmerSeeds(index))?;
        sketch.serialize_field("multi_position_storage", &LegacyMultiStorage(index))?;
        sketch.serialize_field("contigs", &self.contigs)?;
        sketch.serialize_field("total_sequence_length", &self.total_sequence_length)?;
        sketch.serialize_field("contig_lengths", &self.contig_lengths)?;
        sketch.serialize_field("repetitive_kmers", &self.repetitive_kmers)?;
        sketch.serialize_field("marker_seeds", &self.marker_seeds)?;
        sketch.serialize_field("marker_c", &self.marker_c)?;
        sketch.serialize_field("c", &self.c)?;
        sketch.serialize_field("k", &self.k)?;
        sketch.serialize_field("contig_order", &self.contig_order)?;
        sketch.serialize_field("individual_contig", &self.individual_contig)?;
        sketch.serialize_field("amino_acid", &self.amino_acid)?;
        sketch.end()
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
            seed_index: None,
            contigs: vec![],
            total_sequence_length: 0,
            contig_lengths: vec![],
            repetitive_kmers: 0,
            marker_seeds: MMHashSet::default(),
            marker_c: 0,
            c: 0,
            k: 0,
            contig_order:0,
            individual_contig: false,
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
        reverse: bool,
    ) -> Anchor {
        Anchor {
            ref_pos: rpos.0,
            ref_contig: rpos.1,
            query_pos: qpos.0,
            query_contig: qpos.1,
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
