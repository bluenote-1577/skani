use crate::params::*;
//use std::collections::{HashMap, HashSet};
//use std::hash::{BuildHasherDefault, Hash, Hasher};
use smallvec::SmallVec;
use crate::types::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use log::*;

pub fn check_small_contigs(ref_sketches: &Vec<Sketch>, query_sketches: &Vec<Sketch>){
    let mut num_small_sketches = 0;
    let mut num_large_sketches = 0;
    for sketch in ref_sketches.iter(){
        if sketch.marker_seeds.len() < SCREEN_MINIMUM_KMERS{
            num_small_sketches += 1;
        }
        else{
            num_large_sketches += 1;
        }
    }
    for sketch in query_sketches.iter(){
        if sketch.marker_seeds.len() < SCREEN_MINIMUM_KMERS{
            num_small_sketches += 1;
        }
        else{
            num_large_sketches += 1;
        }
    }
    if num_large_sketches + num_small_sketches == 0{
        return
    }
    if num_small_sketches as f64 / (num_large_sketches + num_small_sketches) as f64 > 0.25 &&
    num_small_sketches + num_large_sketches > 10_000{
        log::warn!("Lots of small genomes detected with < 20 marker k-mers. Consider -m or using --faster-small for faster runtimes.");
    }
}

//Used in search, but not in dist,triangle
pub fn screen_refs_indices(
    identity: f64,
    kmer_to_sketch: &KmerToSketch,
    query_sketch: &Sketch,
    sketch_params: &SketchParams,
    ref_sketches: &Vec<Sketch>
) -> Vec<usize>{
    let mut count_hash_map = FxHashMap::default();
    for marker in query_sketch.marker_seeds.iter() {
        if kmer_to_sketch.contains_key(marker) {
            for sketch_id in kmer_to_sketch[marker].iter() {
                let count = count_hash_map.entry(sketch_id).or_insert(0);
                *count += 1;
            }
        }
    }
    //Use fixed K value for AA markers, but flexible ones for DNA because saturation less of an
    //issue.
    let k = if sketch_params.use_aa {
        K_MARKER_AA
    } else {
        K_MARKER_DNA
    };
    let cutoff = identity.powi(k as i32);
    trace!("cutoff screening val {}",cutoff);
    let ret = count_hash_map
        .iter()
        .filter(|x| {
            *x.1 > usize::max((cutoff 
                * usize::min(
                    ref_sketches[**x.0 as usize].marker_seeds.len(),
                    query_sketch.marker_seeds.len(),
                ) as f64) as usize,1)
        })
        .map(|x| **x.0 as usize)
        .collect();
    ret

}

///Quickly check marker k-mers and see
///if the max-contain ANI is > screen_val.
///If rescue_small is enabled, genomes with < 20 markers 
///always pass the filter. Otherwise, at least 1 marker
///needs to be shared between the genomes to return `true`.
pub fn check_markers_quickly(
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
    screen_val: f64,
    rescue_small: bool,
) -> bool{

    if screen_val == 0.{
        return true;
    }
    
    let seeds1;
    let seeds2;
    let min_card;
    if query_sketch.marker_seeds.len() > ref_sketch.marker_seeds.len(){
        seeds1 = &ref_sketch.marker_seeds;
        seeds2 = &query_sketch.marker_seeds;
        min_card = ref_sketch.marker_seeds.len();
    }
    else{
        seeds2 = &ref_sketch.marker_seeds;
        seeds1 = &query_sketch.marker_seeds;
        min_card = query_sketch.marker_seeds.len();
    }
    if min_card < SCREEN_MINIMUM_KMERS && rescue_small{
        return true;
    }

    if min_card == 0{
        if !rescue_small{
            return false;
        }
        else{
            return true;
        }
    }

    assert!(ref_sketch.amino_acid == query_sketch.amino_acid);
    let k = if ref_sketch.amino_acid{K_MARKER_AA} else {K_MARKER_DNA};

    let ratio = screen_val.powi(k.try_into().unwrap()) * min_card as f64;
    let mut ratio = ratio as usize;

    if ratio == 0{
        ratio = 1;
    }

    let mut intersect_len = 0;
    for marker_seed1 in seeds1.iter(){
        if seeds2.contains(marker_seed1){
            intersect_len += 1;
        }
        if intersect_len >= ratio{
            return true;
        }
    }
    trace!("Ratio {}, intersect_len {}, min_card {}", ratio, intersect_len, min_card);
    false
}

///Screen used in triangle, dist, but not search.
///Returns the indices of sketches in ref_sketch that
///pass the filter, using an inverted-index. If `rescue_small` is true,
///query genomes with < 20 k-mers have every reference genome passing the filter.
pub fn screen_refs(
    identity: f64,
    kmer_to_sketch: &KmerToSketch,
    query_sketch: &Sketch,
    sketch_params: &SketchParams,
    ref_sketches: &Vec<Sketch>,
    rescue_small: bool,
) -> FxHashSet<usize> {
    let mut count_hash_map = FxHashMap::default();
    //Don't screen when the input sketch is too small.
    if query_sketch.marker_seeds.len() < 20 && rescue_small{
        return (0..ref_sketches.len()).collect();
    }
    for marker in query_sketch.marker_seeds.iter() {
        if kmer_to_sketch.contains_key(marker) {
            for sketch_id in kmer_to_sketch[marker].iter() {
                let count = count_hash_map.entry(sketch_id).or_insert(0);
                *count += 1;
            }
        }
    }
    //Use fixed K value for AA markers, but flexible ones for DNA because saturation less of an
    //issue.
    let k = if sketch_params.use_aa {
        K_MARKER_AA
    } else {
        K_MARKER_DNA
    };
    let cutoff = identity.powi(k as i32);
    let ret = count_hash_map
        .iter()
        .filter(|x| {
            *x.1 > usize::max((cutoff 
                * usize::min(
                    ref_sketches[**x.0 as usize].marker_seeds.len(),
                    query_sketch.marker_seeds.len(),
                ) as f64) as usize,1)
        })
        .map(|x| **x.0 as usize)
        .collect();
    ret
}
pub fn kmer_to_sketch_from_refs(ref_sketches: &Vec<Sketch>) -> KmerToSketch {
//    let max_size: usize = ref_sketches.iter().map(|x| x.marker_seeds.len()).sum::<usize>();
    let mut ret = KmerToSketch::default();
    //ret.reserve(max_size);
    for (i, ref_sketch) in ref_sketches.iter().enumerate() {
        for kmer in ref_sketch.marker_seeds.iter() {
            let sketch_set = ret.entry(*kmer).or_insert(SmallVec::<[u32; KMER_SK_SMALL_VEC_SIZE]>::new());
            sketch_set.push(i as u32);
        }
    }

//    debug!("{} unique marker k-mers found", ret.len());
//    let mut kmer_stats = vec![];
//    for kmer in ret.keys(){
//        kmer_stats.push(ret[kmer].len());
//    }
//    kmer_stats.sort_unstable();
//    let l = kmer_stats.len();
//    debug!("{} - 10, {} - 50, {} - MAX", kmer_stats[l*1/10], kmer_stats[l/2], kmer_stats[l-1]);
    ret
}
