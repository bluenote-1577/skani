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

    /// Rebuild a finalized index straight from its on-disk parallel arrays.
    ///
    /// The native on-disk layout *is* the in-memory layout, so this is a move of
    /// two already-sorted vectors: no hashing, no sorting, no re-derivation of
    /// the distinct count. This is the whole point of the native format.
    pub fn from_parts(
        seeds: Vec<SeedBits>,
        positions: Vec<SeedPosition>,
        num_distinct_seeds: u32,
    ) -> Self {
        debug_assert_eq!(seeds.len(), positions.len());
        SeedIndex {
            seeds,
            positions,
            build: Vec::new(),
            num_distinct_seeds,
            finalized: true,
        }
    }

    /// The cached distinct-seed count, for the native writer. Unlike
    /// `num_distinct_seeds` this hands back the stored `u32` as-is rather than
    /// widening it, since that is the width it occupies on disk.
    #[inline]
    pub fn distinct_for_io(&self) -> u32 {
        self.num_distinct_seeds
    }

    /// Build from the legacy tagged-hash-map representation found in pre-0.4
    /// `.sketch` files. This is the slow path: a hash map walk followed by a
    /// full sort of every record.
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


//ON-DISK FORMAT VERSIONING
//
//A `.sketch` file, one record inside `sketches.db`, and `markers.bin` are all
//just `bincode::serialize(&(SketchParams, <payload>))`. bincode is not
//self-describing: struct fields are written positionally with no names, no
//tags and no length, so nothing inside a file written by skani <= 0.3.x says
//which layout it is in. Worse, handing such a file to a parser that expects a
//different layout does not fail cleanly - it happily decodes garbage. New
//files therefore carry an explicit header, and old files are recognised by its
//absence.
//
//Layout of a new file: 8 magic bytes, then a u32 little-endian format version,
//then the bincode stream exactly as before. The magic is written as raw bytes
//rather than through bincode so that its on-disk spelling cannot be changed by
//bincode's integer encoding options.
//
//WHY THE MAGIC CANNOT COLLIDE WITH A LEGACY FILE
//
//A legacy stream begins with `SketchParams`, whose first field is `c: usize`.
//bincode's default options encode `usize` as a fixed-width little-endian u64,
//so the first eight bytes of every legacy file are precisely `c` in
//little-endian. That is not a guess: files written by released skani 0.3.1 and
//by 0.3.2 both begin `7d 00 00 00 00 00 00 00`, i.e. c = 125, the default.
//
//`SKETCH_MAGIC` ends in the byte 0xFF, and little-endian puts that last byte in
//the *most significant* byte of the u64. So for a legacy file to open with our
//magic, its `c` would have to be at least 0xFF00_0000_0000_0000, i.e. about
//1.8e19. `c` is the seed compression factor: it is parsed from `-c`, is a few
//hundred in practice, and skani rejects values outside a small range. Even
//ignoring validation entirely, a `c` anywhere near 2^32 would ask for one seed
//every four billion bases, which no genome and no accepted command line can
//produce. The remaining seven magic bytes ("SKANIsk") would have to match on
//top of that. The collision is not improbable, it is unreachable.
pub const SKETCH_MAGIC: [u8; 8] = *b"SKANIsk\xFF";

/// Version of the native on-disk sketch layout this build writes and reads.
pub const SKETCH_FORMAT_VERSION: u32 = 1;

/// Bytes of header in front of the bincode stream of a native-format file.
pub const SKETCH_HEADER_LEN: usize = SKETCH_MAGIC.len() + 4;

/// What the bytes at the front of a serialized sketch stream say about it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SketchFormat {
    /// Current layout. The bincode stream starts at `SKETCH_HEADER_LEN`.
    Native,
    /// No header at all: written by skani <= 0.3.x. The bincode stream starts at
    /// byte 0 and decodes through the legacy on-disk types.
    Legacy,
    /// Our magic, but a version this build does not know. `None` means the
    /// header itself was truncated.
    Unsupported(Option<u32>),
}

/// The native-format header, ready to be written in front of a bincode stream.
pub fn sketch_header_bytes() -> [u8; SKETCH_HEADER_LEN] {
    let mut header = [0u8; SKETCH_HEADER_LEN];
    header[..SKETCH_MAGIC.len()].copy_from_slice(&SKETCH_MAGIC);
    header[SKETCH_MAGIC.len()..].copy_from_slice(&SKETCH_FORMAT_VERSION.to_le_bytes());
    header
}

/// Write the native-format header. Every writer calls this before its bincode
/// stream; nothing writes the legacy layout any more.
pub fn write_sketch_header<W: std::io::Write>(writer: &mut W) -> std::io::Result<()> {
    writer.write_all(&sketch_header_bytes())
}

/// Classify a stream from the bytes at its front. `prefix` may be the whole
/// file, an mmap slice, or just the first `SKETCH_HEADER_LEN` bytes.
pub fn detect_sketch_format(prefix: &[u8]) -> SketchFormat {
    if prefix.len() < SKETCH_MAGIC.len() || prefix[..SKETCH_MAGIC.len()] != SKETCH_MAGIC {
        //No magic. Per the argument above this can only be a pre-0.4 file (or
        //not a sketch at all, in which case the legacy parse fails and the
        //caller reports a corrupt file, exactly as before).
        return SketchFormat::Legacy;
    }
    if prefix.len() < SKETCH_HEADER_LEN {
        return SketchFormat::Unsupported(None);
    }
    let mut version_bytes = [0u8; 4];
    version_bytes.copy_from_slice(&prefix[SKETCH_MAGIC.len()..SKETCH_HEADER_LEN]);
    let version = u32::from_le_bytes(version_bytes);
    if version == SKETCH_FORMAT_VERSION {
        SketchFormat::Native
    } else {
        SketchFormat::Unsupported(Some(version))
    }
}

/// Pull the header off the front of a reader, returning it with the number of
/// bytes actually read. A short read only happens for a file smaller than the
/// header. The bytes are *consumed*, so a caller that then wants to parse the
/// stream as legacy must chain them back on in front of the reader.
pub fn read_sketch_header<R: std::io::Read>(
    reader: &mut R,
) -> std::io::Result<([u8; SKETCH_HEADER_LEN], usize)> {
    let mut header = [0u8; SKETCH_HEADER_LEN];
    let mut filled = 0;
    while filled < header.len() {
        match reader.read(&mut header[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok((header, filled))
}

static LEGACY_FORMAT_WARNED: std::sync::Once = std::sync::Once::new();

/// Emit the pre-0.4 format notice at most once per process.
///
/// The sketch loading paths are rayon-parallel and a user may be loading
/// hundreds of thousands of old sketches, so this has to be both thread-safe
/// and once per *run* rather than once per sketch. `Once` gives us both.
pub fn warn_legacy_sketch_format() {
    LEGACY_FORMAT_WARNED.call_once(|| {
        log::warn!(
            "Loading sketches in the pre-0.4 skani format. They still work and give \
             identical results, but they have to be converted on every load. Re-running \
             `skani sketch` on these genomes rewrites them in the current format and \
             makes loading substantially faster."
        );
    });
}

#[derive(Eq, PartialEq, Deserialize, Clone)]
#[serde(from = "SketchNative")]
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

/// On-disk image of a `SeedIndex` in the native format: the two parallel arrays
/// exactly as they live in memory, plus the cached distinct-seed count. Reading
/// one is two length-prefixed bulk copies, with no hashing and no sort.
#[derive(Deserialize)]
struct SeedIndexNative {
    seeds: Vec<SeedBits>,
    positions: Vec<SeedPosition>,
    num_distinct_seeds: u32,
}

/// Borrowed counterpart of `SeedIndexNative`, so writing copies nothing. A
/// `&[T]` and a `Vec<T>` have the same bincode encoding, so the two agree.
struct SeedIndexNativeRef<'a>(&'a SeedIndex);

impl Serialize for SeedIndexNativeRef<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut index = serializer.serialize_struct("SeedIndexNative", 3)?;
        index.serialize_field("seeds", self.0.seeds())?;
        index.serialize_field("positions", self.0.positions())?;
        index.serialize_field("num_distinct_seeds", &self.0.distinct_for_io())?;
        index.end()
    }
}

/// Writes the `seed_index` field as a real `Option`.
///
/// A marker-only sketch has no seed index at all, and that has to stay
/// distinguishable on disk from a sketch whose index is genuinely empty -
/// inferring it from "seeds, positions and the distinct count are all empty"
/// conflates the two. bincode writes an `Option` as a one-byte tag, so being
/// explicit costs one byte.
struct NativeSeedIndexField<'a>(Option<&'a SeedIndex>);

impl Serialize for NativeSeedIndexField<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self.0 {
            None => serializer.serialize_none(),
            Some(index) => serializer.serialize_some(&SeedIndexNativeRef(index)),
        }
    }
}

/// Native on-disk image of a `Sketch`, used for reading. Mirrors what
/// `impl Serialize for Sketch` writes, field for field.
#[derive(Deserialize)]
struct SketchNative {
    file_name: String,
    seed_index: Option<SeedIndexNative>,
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

impl From<SketchNative> for Sketch {
    fn from(disk: SketchNative) -> Sketch {
        Sketch {
            file_name: disk.file_name,
            seed_index: disk.seed_index.map(|index| {
                SeedIndex::from_parts(index.seeds, index.positions, index.num_distinct_seeds)
            }),
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

/// A `Sketch` read out of the pre-0.4 layout.
///
/// A type can only have one `Deserialize`, and `Sketch`'s is the native one, so
/// the legacy path needs its own handle. Readers pick between
/// `(SketchParams, Sketch)` and `(SketchParams, LegacySketch)` after inspecting
/// the file header; nothing writes this layout any more.
#[derive(Deserialize)]
#[serde(from = "SketchSerde")]
pub struct LegacySketch(pub Sketch);

impl From<SketchSerde> for LegacySketch {
    fn from(disk: SketchSerde) -> LegacySketch {
        LegacySketch(Sketch::from(disk))
    }
}

/// Legacy on-disk image of a `Sketch`, used for reading only.
///
/// This mirrors the pre-flat-index `Sketch` layout field for field, in the same
/// order, with the same types, so deserializing it accepts exactly what skani
/// <= 0.3.2 wrote. The flat `SeedIndex` is then built from the tagged hash map
/// as the last step of the load, which is the cost the native format removes.
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

/// Writes the native layout: the flat seed index as-is, so that loading it is
/// a pair of bulk copies rather than a hash-map build followed by a sort.
///
/// Callers prepend `sketch_header_bytes()`; this impl emits only the bincode
/// body. The written bytes are deterministic, since the index is already in
/// sorted order (the old `HashMap` iteration order was not).
impl Serialize for Sketch {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut sketch = serializer.serialize_struct("Sketch", 13)?;
        sketch.serialize_field("file_name", &self.file_name)?;
        sketch.serialize_field(
            "seed_index",
            &NativeSeedIndexField(self.seed_index.as_ref()),
        )?;
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

/// Result of the chaining DP for one chunk.
///
/// `pointer_vec[i]` is `i`'s best predecessor, or `i` itself when the chain
/// starts there. Because a predecessor is always at a strictly smaller index,
/// the pointers form a forest in which every tree's root is that tree's
/// minimum index - so component membership is recoverable with a single
/// forward sweep and needs no union-find. `comp_vec[i]` is that root;
/// `comp_size` and `comp_best` are indexed *by root* and are meaningless at
/// non-root positions.
pub struct ChainingResult {
    pub pointer_vec: Vec<usize>,
    pub score_vec: Vec<i32>,
    pub comp_vec: Vec<u32>,
    pub comp_size: Vec<u32>,
    pub comp_best: Vec<u32>,
    pub num_chunks: usize,
}

pub struct ChainingResultANI {
    pub pointer_vec: Vec<usize>,
}

#[derive(Eq, PartialEq, PartialOrd, Ord, Debug, Default, Clone, Copy)]
pub struct Anchor {
    pub query_contig: ContigIndex,
    pub query_pos: GnPosition,
    pub ref_contig: ContigIndex,
    pub ref_pos: GnPosition,
    pub reverse_match: bool,
}

#[derive(PartialEq, PartialOrd, Eq, Ord, Debug, Clone, Default)]
pub struct ChainInterval {
    pub score: i32,
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

/// Anchors and query-seed positions split into fragment-sized chunks, held in
/// compressed-sparse-row form: one flat buffer each, plus offsets.
///
/// This used to be `Vec<Vec<Anchor>>` and `Vec<Vec<GnPosition>>`, which cost
/// two heap allocations per chunk and ~250 chunks per comparison. Those ~500
/// allocations per comparison were the whole reason the chunking stage carried
/// an allocator tax under many threads. The anchors were already being copied
/// into the per-chunk vectors, so the flat form copies exactly as much and
/// allocates twice instead of 500 times.
#[derive(Default)]
pub struct AnchorChunks {
    pub anchors: Vec<Anchor>,
    pub chunk_offsets: Vec<u32>,
    pub lengths: Vec<u32>,
    pub seeds: Vec<GnPosition>,
    pub seed_offsets: Vec<u32>,
}

impl AnchorChunks {
    #[inline]
    pub fn num_chunks(&self) -> usize {
        self.chunk_offsets.len().saturating_sub(1)
    }

    #[inline]
    pub fn chunk(&self, i: usize) -> &[Anchor] {
        &self.anchors[self.chunk_offsets[i] as usize..self.chunk_offsets[i + 1] as usize]
    }

    #[inline]
    pub fn seeds_in_chunk(&self, i: usize) -> &[GnPosition] {
        &self.seeds[self.seed_offsets[i] as usize..self.seed_offsets[i + 1] as usize]
    }
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

#[cfg(test)]
mod sketch_format_tests {
    use super::*;

    /// Serialize a sketch the way every writer does: header, then bincode body.
    fn write_native(sketch: &Sketch, params: &SketchParams) -> Vec<u8> {
        let mut bytes = sketch_header_bytes().to_vec();
        bincode::serialize_into(&mut bytes, &(params, sketch)).unwrap();
        bytes
    }

    fn read_native(bytes: &[u8]) -> (SketchParams, Sketch) {
        assert_eq!(detect_sketch_format(bytes), SketchFormat::Native);
        bincode::deserialize(&bytes[SKETCH_HEADER_LEN..]).unwrap()
    }

    fn demo_params() -> SketchParams {
        SketchParams::new(1000, 125, 15, false, false)
    }

    #[test]
    fn header_is_written_and_recognised() {
        let header = sketch_header_bytes();
        assert_eq!(&header[..8], &SKETCH_MAGIC);
        assert_eq!(&header[8..], &SKETCH_FORMAT_VERSION.to_le_bytes());
        assert_eq!(detect_sketch_format(&header), SketchFormat::Native);
    }

    /// The whole safety argument for the magic: the first eight bytes of a
    /// pre-0.4 file are `SketchParams.c` as a little-endian u64, so every value
    /// of `c` that could ever be written must classify as Legacy, never Native.
    #[test]
    fn no_plausible_legacy_c_can_look_native() {
        for c in [
            0u64, 1, 15, 30, 70, 125, 200, 1000, 100_000, u32::MAX as u64,
            // far beyond anything skani accepts, and still not the magic
            1u64 << 40,
            1u64 << 47,
        ] {
            let mut prefix = c.to_le_bytes().to_vec();
            prefix.extend_from_slice(&[0u8; 8]);
            assert_eq!(
                detect_sketch_format(&prefix),
                SketchFormat::Legacy,
                "c = {} was misclassified as a native header",
                c
            );
        }
        // Reading the magic back as a legacy `c` shows how far out of range it is.
        let as_c = u64::from_le_bytes(SKETCH_MAGIC);
        assert!(as_c >= 0xFF00_0000_0000_0000);
    }

    #[test]
    fn truncated_and_future_headers_are_not_mistaken_for_either_format() {
        // Magic present, version field cut off.
        assert_eq!(
            detect_sketch_format(&SKETCH_MAGIC),
            SketchFormat::Unsupported(None)
        );
        // Magic present, version from the future.
        let mut future = SKETCH_MAGIC.to_vec();
        future.extend_from_slice(&(SKETCH_FORMAT_VERSION + 1).to_le_bytes());
        assert_eq!(
            detect_sketch_format(&future),
            SketchFormat::Unsupported(Some(SKETCH_FORMAT_VERSION + 1))
        );
    }

    #[test]
    fn native_round_trip_preserves_the_seed_index() {
        let mut sketch = Sketch::new(1000, 125, 15, "genome.fa".to_string(), false);
        sketch.seed_index = Some(SeedIndex::default());
        sketch.contigs = vec!["ctg1".to_string(), "ctg2".to_string()];
        sketch.contig_lengths = vec![1234, 5678];
        sketch.total_sequence_length = 6912;
        sketch.marker_seeds.insert(42);
        // Pushed out of order and with a repeated seed, so finalize has work to do.
        sketch.add_seed_position(7, SeedPosition::new(900, 1, true));
        sketch.add_seed_position(3, SeedPosition::new(100, 0, false));
        sketch.add_seed_position(7, SeedPosition::new(50, 0, true));
        sketch.add_seed_position(3, SeedPosition::new(10, 0, true));
        sketch.finalize_seed_index();

        let bytes = write_native(&sketch, &demo_params());
        let (_, back) = read_native(&bytes);

        assert!(back == sketch, "native round trip changed the sketch");
        let index = back.seed_index.as_ref().unwrap();
        assert_eq!(index.len(), 4);
        assert_eq!(index.num_distinct_seeds(), 2);
        // The on-disk arrays really are the sorted in-memory ones.
        assert_eq!(index.seeds(), &[3, 3, 7, 7]);
        assert!(index.seeds().windows(2).all(|w| w[0] <= w[1]));
    }

    /// The prototype inferred `seed_index: None` from "seeds, positions and the
    /// distinct count are all empty", which silently turns a genuine zero-seed
    /// sketch into a marker-only one. The explicit Option must keep them apart.
    #[test]
    fn marker_only_and_empty_index_stay_distinct() {
        let params = demo_params();

        let mut marker_only = Sketch::new(1000, 125, 15, "m.fa".to_string(), false);
        marker_only.seed_index = None;
        marker_only.marker_seeds.insert(99);

        let mut empty_index = marker_only.clone();
        empty_index.seed_index = Some(SeedIndex::default());
        empty_index.finalize_seed_index();

        let (_, marker_back) = read_native(&write_native(&marker_only, &params));
        let (_, empty_back) = read_native(&write_native(&empty_index, &params));

        assert!(marker_back.seed_index.is_none(), "marker-only sketch grew an index");
        assert!(empty_back.seed_index.is_some(), "empty index collapsed to None");
        assert_eq!(empty_back.seed_index.as_ref().unwrap().len(), 0);
        assert!(marker_back != empty_back, "the two collapsed into the same sketch");
        // That distinction costs exactly the Option tag byte.
        assert!(marker_back == marker_only);
        assert!(empty_back == empty_index);
    }

    /// A legacy body behind a native header must fail loudly rather than decode
    /// into a plausible-looking sketch.
    #[test]
    fn legacy_bytes_are_not_silently_decoded_as_native() {
        // A legacy stream opens with SketchParams followed by the old Sketch
        // layout, whose second field is `Option<KmerSeeds>` rather than
        // `Option<SeedIndexNative>`. Stand in for it with the params plus a
        // short run of bytes that is not a valid native Sketch.
        let mut body = bincode::serialize(&demo_params()).unwrap();
        body.extend_from_slice(&bincode::serialize(&"genome.fa".to_string()).unwrap());
        body.extend_from_slice(&[1u8]); // Some(..)
        body.extend_from_slice(&u64::MAX.to_le_bytes()); // absurd sequence length
        let parsed: Result<(SketchParams, Sketch), _> = bincode::deserialize(&body);
        assert!(parsed.is_err(), "a legacy-shaped body decoded as native");
    }
}
