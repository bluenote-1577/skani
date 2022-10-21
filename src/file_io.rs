use crate::params::*;
use crate::seeding;
use crate::types::*;
use log::trace;
use log::warn;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

pub fn fastx_to_sketches(ref_files: Vec<String>, sketch_params: &SketchParams) -> Vec<Sketch> {
    let ref_sketches: Mutex<Vec<_>> = Mutex::new(vec![]);
    (0..ref_files.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
            let ref_file = &ref_files[i];
            let mut new_sketch = Sketch::default();
            //    new_sketch.file_name = ref_file.split('/').last().unwrap().to_string();
            new_sketch.file_name = ref_file.to_string();
            new_sketch.c = sketch_params.cs[0];
            let mut i = 0;
            //            let mut reader = parse_fastx_file(&ref_file).expect(ref_file);
            let reader = parse_fastx_file(&ref_file);
            if !reader.is_ok() {
                warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
            } else {
                let mut reader = reader.unwrap();
                trace!("Sketching {}", new_sketch.file_name);
                while let Some(record) = reader.next() {
                    let record = record.expect(&format!("Invalid record for file {}", ref_file));
                    let contig = record.id();
                    let seq = record.seq();
                    if seq.len() >= MIN_LENGTH_CONTIG {
                        new_sketch
                            .contigs
                            .push(String::from_utf8(contig.to_vec()).unwrap());

                        new_sketch.total_sequence_length += seq.len();
                        if sketch_params.use_aa {
                            let orfs = seeding::get_orfs(&seq, sketch_params);
                            if sketch_params.use_syncs {
                                seeding::os_seeds_aa_with_orf(
                                    //            seeding::fmh_seeds_aa(
                                    &seq,
                                    sketch_params,
                                    i as u32,
                                    &mut new_sketch.kmer_seeds_k,
                                    orfs,
                                );
                            } else {
                                seeding::fmh_seeds_aa_with_orf(
                                    //            seeding::fmh_seeds_aa(
                                    &seq,
                                    sketch_params,
                                    i as u32,
                                    &mut new_sketch,
                                    orfs,
                                )
                            }
                        } else {
                            seeding::fmh_seeds(&seq, &sketch_params, i as u32, &mut new_sketch);
                        }
                        i += 1;
                    }
                }
                seeding::get_repetitive_kmers(
                    &new_sketch.kmer_seeds_k,
                    &mut new_sketch.repetitive_kmers,
                );

                {
                    let mut locked = ref_sketches.lock().unwrap();
                    locked.push(new_sketch);
                }
            }
        });
    let mut ref_sketches = ref_sketches.into_inner().unwrap();
    ref_sketches.sort_by(|x, y| x.file_name.cmp(&y.file_name));
    return ref_sketches;
}

pub fn fastx_to_multiple_sketch_rewrite(
    ref_file: &str,
    sketch_params: &SketchParams,
) -> Vec<Sketch> {
    let mut new_sketches = vec![];
    //    let file_name = ref_file.split('/').last().unwrap().to_string();
    let file_name = ref_file.to_string();
    let mut reader = parse_fastx_file(&ref_file).expect("valid path/file");
    println!("Sketching {}", file_name);
    while let Some(record) = reader.next() {
        let mut new_sketch = Sketch::default();
        new_sketch.c = sketch_params.cs[0];
        new_sketch.file_name = file_name.clone();
        let record = record.expect("invalid record");

        let contig = record.id();
        new_sketch
            .contigs
            .push(String::from_utf8(contig.to_vec()).unwrap());
        let seq = record.seq();
        new_sketch.total_sequence_length += seq.len();
        seeding::fmh_seeds(&seq, &sketch_params, 0 as u32, &mut new_sketch);
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
    let mut ani_file = BufWriter::new(File::create(ani_mat_file).expect(file_name));
    let mut af_file = BufWriter::new(File::create(af_mat_file).unwrap());
    write!(&mut ani_file, "{}\n", sketches.len()).unwrap();
    write!(&mut af_file, "{}\n", sketches.len()).unwrap();
    for i in 0..sketches.len() {
        write!(&mut ani_file, "{}", sketches[i].file_name).unwrap();
        write!(&mut af_file, "{}", sketches[i].file_name).unwrap();
        for j in 0..i {
            if anis[j][i].ani == -1. || anis[j][i].ani.is_nan() {
                write!(&mut ani_file, "\t{:.4}", 0.).unwrap();
            } else {
                write!(&mut ani_file, "\t{:.4}", anis[j][i].ani).unwrap();
            }
            write!(&mut af_file, "\t{:.4}", anis[j][i].align_fraction_query).unwrap();
        }
        write!(&mut ani_file, "\n").unwrap();
        write!(&mut af_file, "\n").unwrap();
    }
}

pub fn write_query_ref_list(anis: &Vec<AniEstResult>, file_name: &str) {
    let out_file = format!("{}", file_name);
    let mut out_file = File::create(out_file).expect(file_name);
    for i in 0..anis.len() {
        if anis[i].ani < 0. || anis[i].ani.is_nan() {
            continue;
        }
        let ani = if anis[i].ani < 0. {
            "NA".to_string()
        } else {
            format!("{}", anis[i].ani)
        };
        write!(
            &mut out_file,
            "{}\t{}\t{}\t{}\t{}\n",
            anis[i].ref_file,
            anis[i].query_file,
            ani,
            anis[i].align_fraction_query,
            anis[i].align_fraction_ref
        )
        .unwrap();
    }
}
