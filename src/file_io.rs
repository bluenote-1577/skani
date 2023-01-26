use crate::params::*;

use crate::seeding;
use crate::types::*;
use fxhash::FxHashMap;
use log::*;
use needletail::parse_fastx_file;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::sync::Mutex;

fn write_header(writer: &mut impl Write, id_str: &str, ci: bool) {
    if !ci{
        writeln!(writer,"Ref_file\tQuery_file\t{}\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name", id_str).unwrap();
    } else {
        writeln!(writer,"Ref_file\tQuery_file\t{}\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\t{}_5_percentile\t{}_95_percentile\tStandard_deviation\tRef_mean_ctg_len\tQuery_mean_ctg_len", id_str, id_str, id_str).unwrap();
    }
}

fn write_ani_res(writer: &mut impl Write, ani_res: &AniEstResult, ci: bool) {
    if !ci{
        writeln!(
            writer,
            "{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{}\t{}",
            ani_res.ref_file,
            ani_res.query_file,
            ani_res.ani * 100.,
            ani_res.align_fraction_ref * 100.,
            ani_res.align_fraction_query * 100.,
            ani_res.ref_contig,
            ani_res.query_contig,
        )
        .unwrap();
    } else {
        writeln!(
            writer,
            "{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:0}\t{:0}\t{:0}\t{:0}\t{:0}\t{:0}",
            ani_res.ref_file,
            ani_res.query_file,
            ani_res.ani * 100.,
            ani_res.align_fraction_ref * 100.,
            ani_res.align_fraction_query * 100.,
            ani_res.ref_contig,
            ani_res.query_contig,
            ani_res.num_contigs_r,
            ani_res.num_contigs_q,
            ani_res.ci_lower * 100.,
            ani_res.ci_upper * 100.,
            ani_res.std * 100.,
            ani_res.quant_90_contig_len_r,
            ani_res.quant_50_contig_len_r,
            ani_res.quant_10_contig_len_r,
            ani_res.quant_90_contig_len_q,
            ani_res.quant_50_contig_len_q,
            ani_res.quant_10_contig_len_q,
        )
        .unwrap();
    }
}

pub fn fastx_to_sketches(
    ref_files: &Vec<String>,
    sketch_params: &SketchParams,
    seed: bool,
) -> Vec<Sketch> {
    let ref_sketches: Mutex<Vec<_>> = Mutex::new(vec![]);
    let mut index_vec = (0..ref_files.len()).collect::<Vec<usize>>();
    index_vec.shuffle(&mut thread_rng());
    index_vec.into_par_iter().for_each(|i| {
        let ref_file = &ref_files[i];
        let mut new_sketch = Sketch::new(
            sketch_params.marker_c,
            sketch_params.c,
            sketch_params.k,
            ref_file.to_string(),
            sketch_params.use_aa,
        );
        let reader = parse_fastx_file(ref_file);
        if reader.is_err() {
            if ref_file.contains(".sketch"){
                warn!("{} is not a valid fasta/fastq file but has the .sketch extension. Not all inputs have .sketch extension, so fasta/fastq is assumed.", ref_file);
            }
            else{
                warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
            }
        } else {
            let mut j = 0;
            let mut is_valid = false;
            let mut reader = reader.unwrap();
            trace!("Sketching {} {}", new_sketch.file_name, i);
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record = record.unwrap_or_else(|_| panic!("Invalid record for file {}", ref_file));
                    let contig = record.id();
                    let seq = record.seq();
                    if seq.len() >= MIN_LENGTH_CONTIG {
                        new_sketch
                            .contigs
                            .push(String::from_utf8(contig.to_vec()).unwrap());
                        new_sketch.contig_lengths.push(seq.len() as GnPosition);

                        new_sketch.total_sequence_length += seq.len();
                        if sketch_params.use_aa {
                            let orfs = seeding::get_orfs(&seq, sketch_params);
                            seeding::fmh_seeds_aa_with_orf(
                                &seq,
                                sketch_params,
                                j as u32,
                                &mut new_sketch,
                                orfs,
                                seed,
                            )
                        } else {
                            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                            if is_x86_feature_detected!("avx2") && true{
                                unsafe {
                                    seeding::avx2_fmh_seeds(
                                        &seq,
                                        sketch_params,
                                        j as u32,
                                        &mut new_sketch,
                                        seed,
                                    );
                                }
                            } else {
                                seeding::fmh_seeds(
                                    &seq,
                                    sketch_params,
                                    j as u32,
                                    &mut new_sketch,
                                    seed,
                                );
                            }
                        }
                        //new_sketch.contig_order = 0;
                        j += 1;
                        is_valid = true;
                    }
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                    is_valid = false;
                    break;
                }
            }
            if is_valid {
                if new_sketch.total_sequence_length > 20_000_000 {
                    new_sketch.repetitive_kmers =
                        seeding::get_repetitive_kmers(&new_sketch.kmer_seeds_k);
                }

                {
                    let mut locked = ref_sketches.lock().unwrap();
                    locked.push(new_sketch);
                }
            }
        }
    });
    let mut ref_sketches = ref_sketches.into_inner().unwrap();
    ref_sketches.sort();
    ref_sketches
}
pub fn fastx_to_multiple_sketch_rewrite(
    ref_files: &Vec<String>,
    sketch_params: &SketchParams,
    seed: bool,
) -> Vec<Sketch> {
    let ref_sketches: Mutex<Vec<_>> = Mutex::new(vec![]);
    let mut index_vec = (0..ref_files.len()).collect::<Vec<usize>>();
    index_vec.shuffle(&mut thread_rng());
    index_vec.into_par_iter().for_each(|i| {
        let ref_file = &ref_files[i];
        let reader = parse_fastx_file(ref_file);
        if reader.is_err() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        } else {
            let mut j = 0;
            let mut reader = reader.unwrap();
            trace!("Sketching {} {}", ref_file, i);
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record =
                        record.unwrap_or_else(|_| panic!("Invalid record for file {}", ref_file));
                    let contig = record.id();
                    let seq = record.seq();
                    if seq.len() >= MIN_LENGTH_CONTIG {
                        let mut new_sketch = Sketch::new(
                            sketch_params.marker_c,
                            sketch_params.c,
                            sketch_params.k,
                            ref_file.to_string(),
                            sketch_params.use_aa,
                        );
                        new_sketch
                            .contigs
                            .push(String::from_utf8(contig.to_vec()).unwrap());
                        new_sketch.contig_lengths.push(seq.len() as GnPosition);

                        new_sketch.total_sequence_length += seq.len();
                        if sketch_params.use_aa {
                            let orfs = seeding::get_orfs(&seq, sketch_params);
                            seeding::fmh_seeds_aa_with_orf(
                                &seq,
                                sketch_params,
                                0_u32,
                                &mut new_sketch,
                                orfs,
                                seed,
                            )
                        } else {
                            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                            if is_x86_feature_detected!("avx2") {
                                unsafe {
                                    seeding::avx2_fmh_seeds(
                                        &seq,
                                        sketch_params,
                                        0_u32,
                                        &mut new_sketch,
                                        seed,
                                    );
                                }
                            } else {
                                seeding::fmh_seeds(
                                    &seq,
                                    sketch_params,
                                    0_u32,
                                    &mut new_sketch,
                                    seed,
                                );
                            }
                        }
                        new_sketch.contig_order = j;
                        let mut locked = ref_sketches.lock().unwrap();
                        locked.push(new_sketch);
                        j += 1;
                    }
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                    break;
                }
            }
        }
    });
    let mut ref_sketches = ref_sketches.into_inner().unwrap();
    ref_sketches.sort();
    ref_sketches
}

pub fn write_phyllip_matrix(
    anis: &FxHashMap<usize, FxHashMap<usize, AniEstResult>>,
    sketches: &Vec<Sketch>,
    file_name: &str,
    use_contig_names: bool,
    full_matrix: bool,
    aai: bool,
) {
    let _id_str = if aai { "AAI" } else { "ANI" };
    if file_name.is_empty() {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        writeln!(&mut handle, "{}", sketches.len()).unwrap();
        for i in 0..sketches.len() {
            let name;
            if use_contig_names {
                name = &sketches[i].contigs[0];
            } else {
                name = &sketches[i].file_name;
            }
            write!(&mut handle, "{}", name).unwrap();
            let end;
            if full_matrix {
                end = sketches.len();
            } else {
                end = i;
            }
            for j in 0..end {
                let x = usize::min(i, j);
                let y = usize::max(i, j);
                if i == j {
                    write!(&mut handle, "\t{:.2}", 100.).unwrap();
                } else if !anis.contains_key(&x) || !anis[&x].contains_key(&y) {
                    write!(&mut handle, "\t{:.2}", 0.).unwrap();
                } else if anis[&x][&y].ani == -1. || anis[&x][&y].ani.is_nan() {
                    write!(&mut handle, "\t{:.2}", 0.).unwrap();
                } else {
                    write!(&mut handle, "\t{:.2}", anis[&x][&y].ani * 100.).unwrap();
                }
            }
            writeln!(&mut handle).unwrap();
        }

        let af_mat_file = "skani_matrix.af".to_string();
        let mut af_file = BufWriter::new(File::create(af_mat_file).unwrap());
        writeln!(&mut af_file, "{}", sketches.len()).unwrap();
        for i in 0..sketches.len() {
            let name;
            if use_contig_names {
                name = &sketches[i].contigs[0];
            } else {
                name = &sketches[i].file_name;
            }
            write!(&mut af_file, "{}", name).unwrap();
            let end;
            if full_matrix {
                end = sketches.len();
            } else {
                end = i;
            }
            for j in 0..end {
                if i == j {
                    write!(&mut af_file, "\t{:.2}", 100.).unwrap();
                    continue
                }
                let x = usize::min(i, j);
                let y = usize::max(i, j);
                if !anis.contains_key(&x) || !anis[&x].contains_key(&y) {
                    write!(&mut af_file, "\t{:.2}", 0.).unwrap();
                } else if anis[&x][&y].ani == -1. || anis[&x][&y].ani.is_nan() {
                    write!(&mut af_file, "\t{:.2}", 0.).unwrap();
                } else {
                    if j > i {
                        write!(
                            &mut af_file,
                            "\t{:.2}",
                            anis[&x][&y].align_fraction_ref * 100.
                        )
                        .unwrap();
                    } else {
                        write!(
                            &mut af_file,
                            "\t{:.2}",
                            anis[&x][&y].align_fraction_query * 100.
                        )
                        .unwrap();
                    }
                }
            }
            writeln!(&mut af_file).unwrap();
        }

        info!("Aligned fraction matrix written to skani_matrix.af");
    } else {
        let ani_mat_file = file_name.to_string();
        let af_mat_file = format!("{}.af", file_name);
        let mut ani_file = BufWriter::new(File::create(ani_mat_file).expect(file_name));
        let mut af_file = BufWriter::new(File::create(af_mat_file).unwrap());
        writeln!(&mut ani_file, "{}", sketches.len()).unwrap();
        writeln!(&mut af_file, "{}", sketches.len()).unwrap();
        for i in 0..sketches.len() {
            let name;
            if use_contig_names {
                name = &sketches[i].contigs[0];
            } else {
                name = &sketches[i].file_name;
            }
            write!(&mut ani_file, "{}", name).unwrap();
            write!(&mut af_file, "{}", name).unwrap();
            let end;
            if full_matrix {
                end = sketches.len();
            } else {
                end = i;
            }
            for j in 0..end {
                if i == j{
                    write!(&mut ani_file, "\t{:.2}", 100.).unwrap();
                    write!(&mut af_file, "\t{:.2}", 100.).unwrap();
                    continue
                }
                let x = usize::min(i, j);
                let y = usize::max(i, j);

                if !anis.contains_key(&x) || !anis[&x].contains_key(&y) {
                    write!(&mut ani_file, "\t{:.2}", 0.).unwrap();
                    write!(&mut af_file, "\t{:.2}", 0.).unwrap();
                } else if anis[&x][&y].ani == -1. || anis[&x][&y].ani.is_nan() {
                    write!(&mut ani_file, "\t{:.2}", 0.).unwrap();
                    write!(&mut af_file, "\t{:.2}", 0.).unwrap();
                } else {
                    write!(&mut ani_file, "\t{:.2}", anis[&x][&y].ani * 100.).unwrap();
                    if j > i {
                        write!(
                            &mut af_file,
                            "\t{:.2}",
                            anis[&x][&y].align_fraction_ref * 100.
                        )
                        .unwrap();
                    } else {
                        write!(
                            &mut af_file,
                            "\t{:.2}",
                            anis[&x][&y].align_fraction_query * 100.
                        )
                        .unwrap();
                    }
                }
            }
            writeln!(&mut ani_file).unwrap();
            writeln!(&mut af_file).unwrap();
        }

        info!(
            "Identity and align fraction matrix written to {} and {}.af",
            file_name, file_name
        );
    }
}

pub fn write_sparse_matrix(
    anis: &FxHashMap<usize, FxHashMap<usize, AniEstResult>>,
    _sketches: &Vec<Sketch>,
    file_name: &str,
    aai: bool,
    est_ci: bool,
) {
    let id_str = if aai { "AAI" } else { "ANI" };
    if file_name.is_empty() {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        write_header(&mut handle, id_str, est_ci);
        //        write!(&mut handle,"Ref_file\tQuery_file\t{}\tAlign_fraction_ref\tAlign_fraction_query\t{}_95_percentile\t{}_5_percentile\tRef_name\tQuery_name\n", id_str, id_str, id_str).unwrap();
        for i in anis.keys() {
            for (j, ani_res) in anis[i].iter() {
                if !(anis[i][j].ani == -1. || anis[i][j].ani.is_nan()) {
                    write_ani_res(&mut handle, ani_res, est_ci);
                }
            }
        }
    } else {
        let ani_mat_file = file_name.to_string();
        let mut ani_file = BufWriter::new(File::create(ani_mat_file).expect(file_name));
        write_header(&mut ani_file, id_str, est_ci);
        for i in anis.keys() {
            for (j, ani_res) in anis[i].iter() {
                if !(anis[i][j].ani == -1. || anis[i][j].ani.is_nan()) {
                    write_ani_res(&mut ani_file, ani_res, est_ci);
                }
            }
        }
    }
}

pub fn write_query_ref_list(
    anis: &Vec<AniEstResult>,
    file_name: &str,
    n: usize,
    aai: bool,
    est_ci: bool,
) {
    let id_str = if aai { "AAI" } else { "ANI" };
    let mut query_file_result_map = FxHashMap::default();
    let out_file = file_name.to_string();

    for i in 0..anis.len() {
        if anis[i].ani < 0. || anis[i].ani.is_nan() {
            continue;
        }
        let _ani = if anis[i].ani < 0. {
            "NA".to_string()
        } else {
            format!("{}", anis[i].ani)
        };
        let results = query_file_result_map
            .entry(&anis[i].query_contig)
            .or_insert(vec![]);
        results.push(&anis[i]);
    }
    let mut sorted_keys = query_file_result_map.keys().collect::<Vec<&&String>>();
    sorted_keys.sort();

    if out_file.is_empty() {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        write_header(&mut handle, id_str, est_ci);
        for key in sorted_keys {
            let mut anis = query_file_result_map[key].clone();

            anis.sort_by(|y, x| x.ani.partial_cmp(&y.ani).unwrap());
            for i in 0..usize::min(n, anis.len()) {
                write_ani_res(&mut handle, anis[i], est_ci);
            }
        }
    } else {
        let mut handle;
        handle = File::create(out_file).expect(file_name);
        write_header(&mut handle, id_str, est_ci);
        for key in sorted_keys {
            let mut anis = query_file_result_map[key].clone();

            anis.sort_by(|y, x| x.ani.partial_cmp(&y.ani).unwrap());
            for i in 0..usize::min(n, anis.len()) {
                write_ani_res(&mut handle, anis[i], est_ci);
            }
        }
    }
}

pub fn sketches_from_sketch(ref_files: &Vec<String>) -> (SketchParams, Vec<Sketch>) {
    let ret_sketch_params: Mutex<SketchParams> = Mutex::new(SketchParams::default());
    let ret_ref_sketches: Mutex<Vec<Sketch>> = Mutex::new(vec![]);

    (0..ref_files.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
            let sketch_file = &ref_files[i];
            if !sketch_file.contains("markers.bin") {
                let reader = BufReader::new(File::open(sketch_file).expect(sketch_file));
                let res: Result<(SketchParams, Sketch), _> = bincode::deserialize_from(reader);
                if res.is_ok() {
                    let (temp_sketch_param, temp_ref_sketch) = res.unwrap();
                    let mut locked = ret_sketch_params.lock().unwrap();
                    *locked = temp_sketch_param;
                    let mut locked = ret_ref_sketches.lock().unwrap();
                    locked.push(temp_ref_sketch);
                } else if sketch_file != "markers.bin" {
                    error!(
                        "{} is not a valid .sketch file or is corrupted.",
                        sketch_file
                    );
                }
            }
        });

    let ret_sketch_params = ret_sketch_params.into_inner().unwrap();
    let mut ret_ref_sketches = ret_ref_sketches.into_inner().unwrap();

    ret_ref_sketches.sort_by(|x, y| x.file_name.cmp(&y.file_name));
    (ret_sketch_params, ret_ref_sketches)
}

pub fn marker_sketches_from_marker_file(marker_file: &str) -> (SketchParams, Vec<Sketch>) {
    let reader = BufReader::new(File::open(marker_file).unwrap());
    let res: Result<(SketchParams, Vec<Sketch>), _> = bincode::deserialize_from(reader);
    if res.is_ok() {
        res.unwrap()
    } else {
        error!("Problem reading {}. Exiting. ", marker_file);
        std::process::exit(1)
    }
}
