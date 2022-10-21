use fxhash::FxHashMap;
use fxhash::FxHashSet;
use crate::types::*;
use crate::params::*;

pub fn screen_refs(identity: f64, kmer_to_sketch: &KmerToSketch, query_sketch: &Sketch, sketch_params: &SketchParams) -> Vec<usize>{

    let mut count_hash_map = FxHashMap::default();
    for marker in query_sketch.marker_seeds.iter(){
        if kmer_to_sketch.contains_key(marker){
            for sketch_id in kmer_to_sketch[marker].iter(){
                let count = count_hash_map.entry(sketch_id).or_insert(0);
                *count += 1;
            }
        }
    }
    //Use fixed K value for AA markers, but flexible ones for DNA because saturation less of an
    //issue.
    let k = if sketch_params.use_aa {K_MARKER_AA} else {sketch_params.ks[0]};
    let cutoff = identity.powi(k as i32) * query_sketch.marker_seeds.len() as f64;
    let ret = count_hash_map.iter().filter(|x| *x.1 > (cutoff as usize)).map(|x| **x.0).collect();
    return ret;
}
pub fn kmer_to_sketch_from_refs(ref_sketches: &Vec<Sketch>) -> KmerToSketch {
    let mut ret = KmerToSketch::default();
    for (i,ref_sketch) in ref_sketches.iter().enumerate(){
        for kmer in ref_sketch.marker_seeds.iter(){
            let sketch_set = ret.entry(*kmer).or_insert(FxHashSet::default());
            sketch_set.insert(i);
        }
    }
    return ret
}
