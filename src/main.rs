use clap::{AppSettings, Arg, ArgGroup, Command, SubCommand};
use log::LevelFilter;
use log::{debug, info, warn};
use rayon::prelude::*;
use simple_logging;
use skani::chain;
use skani::file_io;
use skani::params;
use skani::parse;
use skani::screen;
use skani::types;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::sync::Mutex;
use std::time::Instant;
fn main() {
    let matches = Command::new("skani")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1")
        .about("skani")
        .subcommand(
            SubCommand::with_name(params::SKETCH_STRING)
                .arg(
                    Arg::new("reference")
                        .index(1)
                        .help("reference fasta.")
                        .takes_value(true)
                        .multiple(true)
                        .required(true),
                )
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help("compression factor.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .help("use amino acid alphabet."),
                )
                .arg(Arg::new("o").short('o').help("output. ").takes_value(true))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .arg(Arg::new("v").short('v').help("verbose."))
                .arg(Arg::new("trace").long("trace").help("verbose.")),
        )
        .subcommand(
            SubCommand::with_name(params::DIST_STRING)
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help("compression factor.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("query")
                        .short('q')
                        .help("query fasta(s) or sketches(es).")
                        .takes_value(true)
                        .required(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .help("reference fasta(s) or sketch(es).")
                        .takes_value(true)
                        .required(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .help("use amino acid alphabet."),
                )
                .arg(Arg::new("v").short('v').help("verbose."))
                .arg(Arg::new("trace").long("trace").help("verbose."))
                .arg(Arg::new("o").short('o').help("output. ").takes_value(true))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .arg(Arg::new("s").short('s').takes_value(true).help("screen. "))
                .arg(Arg::new("robust").long("robust").help("robust ani estimation; trim off 5/95% quantiles. "))
        )
        .subcommand(
            SubCommand::with_name(params::TRIANGLE_STRING)
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help("compression factor.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("fasta_files")
                        .index(1)
                        .help("fasta(s) or sketch(es).")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("fasta_list")
                        .short('l')
                        .help("file with each line containing one fasta file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .help("use amino acid alphabet."),
                )
                .arg(Arg::new("v").short('v').help("verbose."))
                .arg(Arg::new("trace").long("trace").help("verbose."))
                .arg(Arg::new("o").short('o').help("output. ").takes_value(true))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .arg(Arg::new("s").short('s').takes_value(true).help("screen. "))
                .group(
                    ArgGroup::new("ref")
                        .arg("fasta_files")
                        .arg("fasta_list")
                        .required(true),
                )
                .arg(Arg::new("robust").long("robust").help("robust ani estimation; trim off 5/95% quantiles. "))
        )
        .get_matches();

    let (mut sketch_params, command_params) = parse::parse_params(&matches);
    if command_params.mode == params::Mode::Sketch {
        let ref_sketches = file_io::fastx_to_sketches(command_params.ref_files, &sketch_params);
        let kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches);
        let mut file_bin = BufWriter::new(File::create(command_params.out_file_name).unwrap());
        bincode::serialize_into(
            &mut file_bin,
            &(sketch_params, ref_sketches, kmer_to_sketch),
        )
        .unwrap();
    } else {
        let now = Instant::now();
        let ref_sketches;
        let kmer_to_sketch;

        if command_params.ref_is_sketch {
            let reader = BufReader::new(File::open(&command_params.ref_files[0]).unwrap());
            let (temp_sketch_param, temp_ref_sketches, temp_kmer_to_sketch): (
                params::SketchParams,
                Vec<types::Sketch>,
                types::KmerToSketch,
            ) = bincode::deserialize_from(reader).unwrap();
            ref_sketches = temp_ref_sketches;
            kmer_to_sketch = temp_kmer_to_sketch;
            if sketch_params != temp_sketch_param {
                warn!("Input sketch parameters different than reference sketch parameters; using reference sketch parameters");
            }
            sketch_params = temp_sketch_param;
        } else {
            ref_sketches = file_io::fastx_to_sketches(command_params.ref_files, &sketch_params);
            kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches)
        }

        let query_sketches;
        if command_params.query_is_sketch {
            let reader = BufReader::new(File::open(&command_params.query_files[0]).unwrap());
            let (temp_query_sketch_param, temp_query_sketches, _): (
                params::SketchParams,
                Vec<types::Sketch>,
                types::KmerToSketch,
            ) = bincode::deserialize_from(reader).unwrap();
            query_sketches = temp_query_sketches;
            if sketch_params != temp_query_sketch_param {
                panic!("Query sketch parameters were not equal to reference sketch parameters. Exiting.");
            }
        } else {
            query_sketches = file_io::fastx_to_sketches(command_params.query_files, &sketch_params);
        }

        info!("Generating sketch time: {}", now.elapsed().as_secs_f32());
        let now = Instant::now();
        if command_params.mode == params::Mode::Dist {
            let js = (0..query_sketches.len()).collect::<Vec<usize>>();
            let anis: Mutex<Vec<_>> = Mutex::new(vec![]);
            js.into_par_iter().for_each(|j| {
                let query_sketch = &query_sketches[j];
                let screened_refs;
                if command_params.screen {
                    screened_refs = screen::screen_refs(
                        command_params.screen_val,
                        &kmer_to_sketch,
                        query_sketch,
                        &sketch_params,
                    );
                    info!(
                        "{} has {} refs passing screening.",
                        query_sketch.file_name,
                        screened_refs.len()
                    );
                } else {
                    screened_refs = (0..ref_sketches.len()).collect::<Vec<usize>>();
                }
                screened_refs.into_par_iter().for_each(|i| {
                    let ref_sketch = &ref_sketches[i];
                    let map_params =
                        chain::map_params_from_sketch(ref_sketch, sketch_params.use_aa, command_params.robust);
                    let ani_res;
                    if map_params != params::MapParams::default() {
                        ani_res = chain::chain_seeds(ref_sketch, &query_sketch, map_params);
                    }
                    else{
                        ani_res = types::AniEstResult::default();
                    }
                    let mut locked = anis.lock().unwrap();
                    locked.push(ani_res);
                });
            });
            let anis = anis.into_inner().unwrap();
            file_io::write_query_ref_list(&anis, &command_params.out_file_name);
        } else if command_params.mode == params::Mode::Triangle {
            let anis: Mutex<Vec<_>> =
                Mutex::new(vec![
                    vec![types::AniEstResult::default(); ref_sketches.len()];
                    ref_sketches.len()
                ]);
            (0..ref_sketches.len() - 1)
                .collect::<Vec<usize>>()
                .into_par_iter()
                .for_each(|i| {
                    let screened_refs;
                    let ref_sketch_i = &ref_sketches[i];
                    if command_params.screen {
                        screened_refs = screen::screen_refs(
                            command_params.screen_val,
                            &kmer_to_sketch,
                            &ref_sketch_i,
                            &sketch_params,
                        );
                        info!(
                            "{} has {} refs passing screening.",
                            ref_sketch_i.file_name,
                            screened_refs.len()
                        );
                    } else {
                        screened_refs = (i + 1..ref_sketches.len()).collect::<Vec<usize>>();
                    }
                    for j in i + 1..ref_sketches.len() {
                        if screened_refs.contains(&j) {
                            let map_params =
                                chain::map_params_from_sketch(ref_sketch_i, sketch_params.use_aa, command_params.robust);
                            let ref_sketch_j = &ref_sketches[j];
                            let ani_res =
                                chain::chain_seeds(ref_sketch_i, ref_sketch_j, map_params);
                            let mut locked = anis.lock().unwrap();
                            locked[i][j] = ani_res;
                        }
                    }
                });
            let anis = anis.into_inner().unwrap();
            file_io::write_phyllip_matrix(&anis, &ref_sketches, &command_params.out_file_name);
        }
        info!("Alignment time: {}", now.elapsed().as_secs_f32());
    }
}
