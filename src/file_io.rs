use crate::params::*;
use crate::seeding;
use crate::types::*;
use needletail::{parse_fastx_file};

pub fn fastx_to_sketch_rewrite(
    ref_file: &str,
    ks: &Vec<usize>,
    cs: &Vec<usize>,
    use_syncs: bool,
) -> Sketch {
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
        let seq_len = seq.len();
//        if seq_len < {
            //                                println!("Sequence {} skipped because the length {} is less than {}", contig, seq.len(), TRIANGLE_LENGTH_CUTOFF);
//            continue;
//        }
        new_sketch.total_sequence_length += seq.len();
        if use_syncs {
//            let adj_c = seeding::open_sync_seeds(&seq, k, c, i as u32, &mut new_sketch.kmer_seeds);
            seeding::fmh_seeds(&seq, ks, cs, i as u32, &mut new_sketch.kmer_seeds_k);
//            new_sketch.c_adj = adj_c;
        } else {
            seeding::fmh_seeds(&seq, ks, cs, i as u32, &mut new_sketch.kmer_seeds_k);
//            new_sketch.c_adj = adj_c;
        }
        i += 1;
    }

    return new_sketch
}

pub fn fastx_to_multiple_sketch_rewrite(
    ref_file: &str,
    k: &Vec<usize>,
    c: &Vec<usize>,
    use_syncs: bool,
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
//        if seq_len < TRIANGLE_LENGTH_CUTOFF {
            //                                println!("Sequence {} skipped because the length {} is less than {}", contig, seq.len(), TRIANGLE_LENGTH_CUTOFF);
//            continue;
//        }
        new_sketch.total_sequence_length += seq.len();
        if use_syncs {
//            let adj_c = seeding::open_sync_seeds(&seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds);
            seeding::fmh_seeds(&seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds_k);
//            new_sketch.c_adj = adj_c;
        } else {
            seeding::fmh_seeds(&seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds_k);
//            new_sketch.c_adj = adj_c;
        }
        new_sketches.push(new_sketch);
    }

    return new_sketches
}

