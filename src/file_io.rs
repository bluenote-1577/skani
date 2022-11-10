use crate::params::*;
use rand::thread_rng;
use rand::seq::SliceRandom;
use crate::seeding;
use crate::types::*;
use bio::io::fasta;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use log::trace;
use log::warn;
use needletail::parse_fastx_file;
use pariter::IteratorExt as _;
use rayon::prelude::*;
use seq_io::fastx::Reader;
use seq_io::BaseRecord;
use smallvec::SmallVec;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::sync::Mutex;

pub fn fastx_to_sketches(
    ref_files: &Vec<String>,
    sketch_params: &SketchParams,
    seed: bool,
) -> Vec<Sketch> {
    let ref_sketches: Mutex<Vec<_>> = Mutex::new(vec![]);

    let mut index_vec = (0..ref_files.len()).collect::<Vec<usize>>();
    index_vec.shuffle(&mut thread_rng());
    index_vec
        .into_par_iter()
        .for_each(|i| {
            let ref_file = &ref_files[i];
            let mut new_sketch = Sketch::default();
            //    new_sketch.file_name = ref_file.split('/').last().unwrap().to_string();
            new_sketch.file_name = ref_file.to_string();
            new_sketch.c = sketch_params.c;
            new_sketch.k = sketch_params.k;
            //            let mut reader = Reader::from_path(ref_file);
            let reader = parse_fastx_file(&ref_file);
            if !reader.is_ok() {
                warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
            } else {
                let mut j = 0;
                let mut is_valid = true;
                let mut reader = reader.unwrap();
                trace!("Sketching {} {}", new_sketch.file_name, i);
                while let Some(record) = reader.next() {
                    if record.is_ok() {
                        let record =
                            record.expect(&format!("Invalid record for file {}", ref_file));
                        let contig = record.id();
                        let seq = record.seq();
                        if seq.len() >= MIN_LENGTH_CONTIG {
                            new_sketch
                                .contigs
                                .push(String::from_utf8(contig.to_vec()).unwrap());
                            new_sketch.contig_lengths.push(seq.len() as GnPosition);
                            //                                .push(contig.unwrap().to_string());

                            new_sketch.total_sequence_length += seq.len();
                            if sketch_params.use_aa {
                                let orfs = seeding::get_orfs(&seq, sketch_params);
                                seeding::fmh_seeds_aa_with_orf(
                                    //            seeding::fmh_seeds_aa(
                                    &seq,
                                    sketch_params,
                                    j as u32,
                                    &mut new_sketch,
                                    orfs,
                                    seed,
                                )
                            } else {
                                seeding::fmh_seeds(
                                    &seq,
                                    &sketch_params,
                                    j as u32,
                                    &mut new_sketch,
                                    seed,
                                );
                            }
                            j += 1;
                        }
                    } else {
                        warn!("File {} is not a valid fasta/fastq file", ref_file);
                        is_valid = false;
                        break;
                    }
                }
                if is_valid {
                    if new_sketch.total_sequence_length > 20000000 {
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
    ref_sketches.sort_by(|x, y| x.file_name.cmp(&y.file_name));
    return ref_sketches;
}
//NOT FUNCTIONAL
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
        new_sketch.c = sketch_params.c;
        new_sketch.file_name = file_name.clone();
        let record = record.expect("invalid record");

        let contig = record.id();
        new_sketch
            .contigs
            .push(String::from_utf8(contig.to_vec()).unwrap());
        let seq = record.seq();
        new_sketch.total_sequence_length += seq.len();
        seeding::fmh_seeds(&seq, &sketch_params, 0 as u32, &mut new_sketch, false);
        new_sketches.push(new_sketch);
    }

    return new_sketches;
}

pub fn write_phyllip_matrix(
    anis: &FxHashMap<usize, FxHashMap<usize, AniEstResult>>,
    sketches: &Vec<Sketch>,
    file_name: &str,
) {
    let file_name_wr;
    if file_name == "" {
        file_name_wr = "skani_res"
    } else {
        file_name_wr = file_name
    }
    let ani_mat_file = format!("{}.ani", file_name_wr);
    let af_mat_file = format!("{}.af", file_name_wr);
    let mut ani_file = BufWriter::new(File::create(ani_mat_file).expect(file_name_wr));
    let mut af_file = BufWriter::new(File::create(af_mat_file).unwrap());
    write!(&mut ani_file, "{}\n", sketches.len()).unwrap();
    write!(&mut af_file, "{}\n", sketches.len()).unwrap();
    for i in 0..sketches.len() {
        write!(&mut ani_file, "{}", sketches[i].file_name).unwrap();
        write!(&mut af_file, "{}", sketches[i].file_name).unwrap();
        for j in 0..i {
            if !anis.contains_key(&j) || !anis[&j].contains_key(&i) {
                write!(&mut ani_file, "\t{:.4}", 0.).unwrap();
                write!(&mut af_file, "\t{:.4}", 0.).unwrap();
            } else if anis[&j][&i].ani == -1. || anis[&j][&i].ani.is_nan() {
                write!(&mut ani_file, "\t{:.4}", 0.).unwrap();
                write!(&mut af_file, "\t{:.4}", 0.).unwrap();
            } else {
                write!(&mut ani_file, "\t{:.4}", anis[&j][&i].ani).unwrap();
                write!(&mut af_file, "\t{:.4}", anis[&j][&i].align_fraction_query).unwrap();
            }
        }
        write!(&mut ani_file, "\n").unwrap();
        write!(&mut af_file, "\n").unwrap();
    }
}

pub fn write_sparse_matrix(
    anis: &FxHashMap<usize, FxHashMap<usize, AniEstResult>>,
    sketches: &Vec<Sketch>,
    file_name: &str,
) {
    let file_name_wr;
    if file_name == "" {
        file_name_wr = "skani_res"
    } else {
        file_name_wr = file_name
    }
    let ani_mat_file = format!("{}.sparse", file_name_wr);
    let mut ani_file = BufWriter::new(File::create(ani_mat_file).expect(file_name_wr));
    for i in anis.keys() {
        for (j, ani_res) in anis[i].iter() {
            if !(anis[i][j].ani == -1. || anis[i][j].ani.is_nan()) {
                write!(
                    &mut ani_file,
                    "{}\t{}",
                    sketches[*i].file_name, sketches[*j].file_name
                )
                .unwrap();
                write!(
                    &mut ani_file,
                    "\t{:.4}\t{:.4}\n",
                    ani_res.ani, ani_res.align_fraction_query
                )
                .unwrap();
            }
        }
    }
}

pub fn write_query_ref_list(anis: &Vec<AniEstResult>, file_name: &str, n: usize) {
    let mut query_file_result_map = FxHashMap::default();
    let out_file = format!("{}", file_name);

    for i in 0..anis.len() {
        if anis[i].ani < 0. || anis[i].ani.is_nan() {
            continue;
        }
        let ani = if anis[i].ani < 0. {
            "NA".to_string()
        } else {
            format!("{}", anis[i].ani)
        };
        let results = query_file_result_map
            .entry(&anis[i].query_file)
            .or_insert(vec![]);
        results.push(&anis[i]);
    }
    let mut sorted_keys = query_file_result_map.keys().collect::<Vec<&&String>>();
    sorted_keys.sort();

    if out_file == "" {
        let stdout = io::stdout();
        let mut handle = stdout.lock();

        write!(&mut handle,"Ref_file\tQuery_file\tANI/AAI\tAlign_fraction_query\tAlign_fraction_reference\tANI/AAI_95_percentile\tANI/AAI_5_percentile\tRef_name\tQuery_name\n").unwrap();
        for key in sorted_keys {
            let mut anis = query_file_result_map[key].clone();

            anis.sort_by(|y, x| x.ani.partial_cmp(&y.ani).unwrap());
            for i in 0..usize::min(n, anis.len()) {
                write!(
                    &mut handle,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    anis[i].ref_file,
                    anis[i].query_file,
                    anis[i].ani,
                    anis[i].align_fraction_query,
                    anis[i].align_fraction_ref,
                    anis[i].ci_upper,
                    anis[i].ci_lower,
                    anis[i].ref_contig,
                    anis[i].query_contig,
                )
                .unwrap();
            }
        }
    } else {
        let mut handle;
        handle = File::create(out_file).expect(file_name);
        write!(&mut handle,"Ref_file\tQuery_file\tANI/AAI\tAlign_fraction_query\tAlign_fraction_reference\tANI/AAI_95_percentile\tANI/AAI_5_percentile\n").unwrap();
        for key in sorted_keys {
            let mut anis = query_file_result_map[key].clone();

            anis.sort_by(|y, x| x.ani.partial_cmp(&y.ani).unwrap());
            for i in 0..usize::min(n, anis.len()) {
                write!(
                    &mut handle,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    anis[i].ref_file,
                    anis[i].query_file,
                    anis[i].ani,
                    anis[i].align_fraction_query,
                    anis[i].align_fraction_ref,
                    anis[i].ci_upper,
                    anis[i].ci_lower,
                )
                .unwrap();
            }
        }
    }
}

pub fn sketches_from_sketch(
    ref_files: &Vec<String>,
    marker: bool,
) -> (SketchParams, Vec<Sketch>) {
    let ret_sketch_params: Mutex<SketchParams> = Mutex::new(SketchParams::default());
    let ret_ref_sketches: Mutex<Vec<Sketch>> = Mutex::new(vec![]);
    let mut test_map = KmerToSketch::default();

    (0..ref_files.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
        let sketch_file = &ref_files[i];
        if marker && sketch_file.contains(".marker") ||
        !marker && !sketch_file.contains(".marker"){
            let reader = BufReader::new(File::open(sketch_file).expect(sketch_file));
            let (temp_sketch_param, temp_ref_sketch): (
                SketchParams,
                Sketch,
            ) = bincode::deserialize_from(reader).unwrap();
            //let temp_kmer_to_sketch = temp_kmer_to_sketch.into_iter().map(|x| (x.0,x.1.into_iter().collect::<FxHashSet<usize>>())).collect::<KmerToSketch>();
            let mut locked = ret_sketch_params.lock().unwrap();
            *locked = temp_sketch_param;
            let mut locked = ret_ref_sketches.lock().unwrap();
            locked.push(temp_ref_sketch);
        }
    });

    let ret_sketch_params = ret_sketch_params.into_inner().unwrap();
    let mut ret_ref_sketches = ret_ref_sketches.into_inner().unwrap();

    ret_ref_sketches.sort_by(|x, y| x.file_name.cmp(&y.file_name));
    return (ret_sketch_params, ret_ref_sketches);
}

pub fn seed_screened_sketches(
    sketch_params: &SketchParams,
    ref_sketches: &mut Vec<Sketch>,
    screened_refs: &FxHashSet<usize>,
) {
    let file_names = screened_refs
        .iter()
        .map(|x| (*x, ref_sketches[*x].file_name.clone()))
        .collect::<FxHashMap<usize, String>>();
    let mutex_refs: Mutex<&mut Vec<Sketch>> = Mutex::new(ref_sketches);
    screened_refs.into_par_iter().for_each(|i| {
        let ref_file = &file_names[i];
        let mut new_sketch = Sketch::default();
        new_sketch.file_name = ref_file.to_string();
        new_sketch.c = sketch_params.c;
        new_sketch.k = sketch_params.k;

        let reader = parse_fastx_file(&ref_file);
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        } else {
            let mut is_valid = true;
            let mut reader = reader.unwrap();
            let mut j = 0;
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record = record.expect(&format!("Invalid record for file {}", ref_file));
                    let contig = record.id();
                    let seq = record.seq();
                    if seq.len() >= MIN_LENGTH_CONTIG {
                        new_sketch.total_sequence_length += seq.len();
                        if sketch_params.use_aa {
                            let orfs = seeding::get_orfs(&seq, sketch_params);
                            seeding::fmh_seeds_aa_with_orf(
                                //            seeding::fmh_seeds_aa(
                                &seq,
                                sketch_params,
                                j as u32,
                                &mut new_sketch,
                                orfs,
                                true,
                            )
                        } else {
                            seeding::fmh_seeds(
                                &seq,
                                &sketch_params,
                                j as u32,
                                &mut new_sketch,
                                true,
                            );
                        }
                        j += 1;
                    }
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                    is_valid = false;
                    break;
                }
            }
            let mut locked = mutex_refs.lock().unwrap();
            locked[*i] = new_sketch;
        }
    });
}
