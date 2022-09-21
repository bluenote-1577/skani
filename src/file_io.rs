use crate::params::*;
use crate::seeding;
use crate::types::*;
use needletail::parse_fastx_file;

pub fn fastx_to_sketch_rewrite(ref_file: &str, sketch_params: &SketchParams) -> Sketch {
    let mut new_sketch = Sketch::default();
    new_sketch.file_name = ref_file.split('/').last().unwrap().to_string();
    let mut i = 0;
    let mut reader = parse_fastx_file(&ref_file).expect("valid path/file");
    while let Some(record) = reader.next() {
        let record = record.expect("invalid record");

        let contig = record.id();
        new_sketch
            .contigs
            .push(String::from_utf8(contig.to_vec()).unwrap());
        let seq = record.seq();
        new_sketch.total_sequence_length += seq.len();
        if sketch_params.use_aa {
            seeding::fmh_seeds_aa(
                &seq,
                &sketch_params.ks,
                &sketch_params.cs,
                i as u32,
                &mut new_sketch.kmer_seeds_k,
            );
        } else {
            seeding::fmh_seeds(
                &seq,
                &sketch_params.ks,
                &sketch_params.cs,
                i as u32,
                &mut new_sketch.kmer_seeds_k,
            );
        }
        i += 1;
    }

    return new_sketch;
}

pub fn fastx_to_multiple_sketch_rewrite(
    ref_file: &str,
    sketch_params: &SketchParams,
) -> Vec<Sketch> {
    let mut new_sketches = vec![];
    let file_name = ref_file.split('/').last().unwrap().to_string();
    let mut reader = parse_fastx_file(&ref_file).expect("valid path/file");
    while let Some(record) = reader.next() {
        let mut new_sketch = Sketch::default();
        new_sketch.file_name = file_name.clone();
        let record = record.expect("invalid record");

        let contig = record.id();
        new_sketch
            .contigs
            .push(String::from_utf8(contig.to_vec()).unwrap());
        let seq = record.seq();
        new_sketch.total_sequence_length += seq.len();
        seeding::fmh_seeds(
            &seq,
            &sketch_params.ks,
            &sketch_params.cs,
            0 as u32,
            &mut new_sketch.kmer_seeds_k,
        );
        new_sketches.push(new_sketch);
    }

    return new_sketches;
}
