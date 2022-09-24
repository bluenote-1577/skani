use crate::params::*;
use crate::seeding;
use crate::types::*;
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::Write;

pub fn fastx_to_sketch_rewrite(ref_file: &str, sketch_params: &SketchParams) -> Sketch {
    let mut new_sketch = Sketch::default();
    new_sketch.file_name = ref_file.split('/').last().unwrap().to_string();
    let mut i = 0;
    let mut reader = parse_fastx_file(&ref_file).expect(ref_file);
    println!("Sketching {}", new_sketch.file_name);
    while let Some(record) = reader.next() {
        let record = record.expect("invalid record");

        let contig = record.id();
        new_sketch
            .contigs
            .push(String::from_utf8(contig.to_vec()).unwrap());
        let seq = record.seq();
        new_sketch.total_sequence_length += seq.len();
        if sketch_params.use_aa {
            let orfs = seeding::get_orfs(&seq, sketch_params);
            seeding::fmh_seeds_aa_with_orf(
                //            seeding::fmh_seeds_aa(
                &seq,
                sketch_params,
                i as u32,
                &mut new_sketch.kmer_seeds_k,
                orfs,
            );
            seeding::get_repetitive_kmers(
                &new_sketch.kmer_seeds_k,
                &mut new_sketch.repetitive_kmers,
            );
        } else {
            seeding::fmh_seeds(
                &seq,
                &sketch_params.ks,
                &sketch_params.cs,
                i as u32,
                &mut new_sketch.kmer_seeds_k,
            );
            seeding::get_repetitive_kmers(
                &new_sketch.kmer_seeds_k,
                &mut new_sketch.repetitive_kmers,
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
    println!("Sketching {}", file_name);
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

pub fn write_phyllip_matrix(
    anis: &Vec<Vec<AniEstResult>>,
    sketches: &Vec<Sketch>,
    file_name: &str,
) {
    let ani_mat_file = format!("{}.ani", file_name);
    let af_mat_file = format!("{}.af", file_name);
    let mut ani_file = File::create(ani_mat_file).expect(file_name);
    let mut af_file = File::create(af_mat_file).unwrap();
    write!(&mut ani_file, "{}\n", sketches.len()).unwrap();
    write!(&mut af_file, "{}\n", sketches.len()).unwrap();
    for i in 0..sketches.len() {
        write!(&mut ani_file, "{}", sketches[i].file_name).unwrap();
        write!(&mut af_file, "{}", sketches[i].file_name).unwrap();
        for j in 0..i {
            write!(&mut ani_file, "\t{:.4}", anis[j][i].ani).unwrap();
            write!(&mut af_file, "\t{:.4}", anis[j][i].align_fraction).unwrap();
        }
        write!(&mut ani_file, "\n").unwrap();
        write!(&mut af_file, "\n").unwrap();
    }
}
