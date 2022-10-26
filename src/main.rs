use clap::{AppSettings, Arg, ArgGroup, Command, SubCommand};
use fxhash::FxHashMap;
use fxhash::FxHashSet;
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
                    Arg::new("fasta_files")
                        .index(1)
                        .help("reference fasta.")
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
                .group(
                    ArgGroup::new("ref")
                        .arg("fasta_files")
                        .arg("fasta_list")
                        .required(true),
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
                .arg(
                    Arg::new("s")
                        .short('s')
                        .long("screen")
                        .help("sketch for screen; don't seed."),
                )
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
                        .multiple(true),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .help("reference fasta(s) or sketch(es).")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("ref_fasta_list")
                        .long("rl")
                        .help("file with each line containing one fasta file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("query_fasta_list")
                        .long("ql")
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
                .group(
                    ArgGroup::new("ref")
                        .arg("reference")
                        .arg("ref_fasta_list")
                        .required(true),
                )
                .group(
                    ArgGroup::new("q")
                        .arg("query")
                        .arg("query_fasta_list")
                        .required(true),
                )
                .arg(Arg::new("s").short('s').takes_value(true).help("screen. "))
                .arg(
                    Arg::new("robust")
                        .long("robust")
                        .help("robust ani estimation; trim off 5/95% quantiles. "),
                )
                .arg(
                    Arg::new("median")
                        .long("median")
                        .help("median ani estimation."),
                ),
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
                .arg(
                    Arg::new("robust")
                        .long("robust")
                        .help("robust ani estimation; trim off 5/95% quantiles. "),
                )
                .arg(
                    Arg::new("median")
                        .long("median")
                        .help("median ani estimation."),
                )
                .arg(
                    Arg::new("sparse")
                        .long("sparse")
                        .help("output sparse matrix."),
                ),
        )
        .get_matches();

    let (mut sketch_params, command_params) = parse::parse_params(&matches);
    if command_params.ref_files.len() == 0 {
        panic!("No reference fastas/sketches found.")
    }
    let all_range_as_set = (0..command_params.ref_files.len()).collect::<FxHashSet<usize>>();

    //SKETCHING
    if command_params.mode == params::Mode::Sketch {
        let ref_sketches = file_io::fastx_to_sketches(
            &command_params.ref_files,
            &sketch_params,
            !command_params.screen,
        );
        //        let kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches);
        //        let mut flat_kmer_to_sketch: Vec<(types::KmerBits, Vec<usize>)> =
        //            Vec::with_capacity(kmer_to_sketch.len());
        //        for (key, set) in kmer_to_sketch {
        //            let vec = set.iter().map(|x| *x).collect::<Vec<usize>>();
        //            flat_kmer_to_sketch.push((key, vec));
        //        }
        let len = 100;
        let num_iters = ref_sketches.len() / len + 1;
        for i in 0..num_iters {
            let mut file_bin = BufWriter::new(
                File::create(format!("{}{}", &command_params.out_file_name, i)).unwrap(),
            );
            let ind1 = i * len;
            let ind2 = usize::min(ref_sketches.len(), (i + 1) * len);

            //            let ind1_flat = i * flat_kmer_to_sketch.len()/num_iters;
            //            let ind2_flat = usize::min(
            //                flat_kmer_to_sketch.len(),
            //                (i + 1) * flat_kmer_to_sketch.len() / num_iters,
            //            );
            bincode::serialize_into(
                &mut file_bin,
                &(
                    &sketch_params,
                    &ref_sketches[ind1..ind2],
                    //                    &flat_kmer_to_sketch[ind1_flat..ind2_flat],
                ),
            )
            .unwrap();
        }
        return;
    }
    //

    let now = Instant::now();
    let ref_sketches;
    let kmer_to_sketch;

    if command_params.refs_are_sketch {
        info!("Sketches detected.");
        (sketch_params, ref_sketches) =
            file_io::sketches_from_sketch(&command_params.ref_files, &sketch_params);
        if command_params.screen {
            kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches);
        } else {
            kmer_to_sketch = types::KmerToSketch::default();
        }

        if sketch_params.seed == false {
            if command_params.mode != params::Mode::Dist || !command_params.screen {
                panic!("Sketched references are not seeded; must use `skani dist` with screening parameter `-s`.");
            }
        }
    } else {
        if command_params.screen && command_params.mode == params::Mode::Dist {
            ref_sketches =
                file_io::fastx_to_sketches(&command_params.ref_files, &sketch_params, false);
        } else {
            ref_sketches =
                file_io::fastx_to_sketches(&command_params.ref_files, &sketch_params, true);
        }
        if command_params.screen{
            kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches)
        }
        else{
            kmer_to_sketch = types::KmerToSketch::default();
        }
    }

    let query_sketches;
    let query_params;
    if command_params.queries_are_sketch {
        (query_params, query_sketches) =
            file_io::sketches_from_sketch(&command_params.query_files, &sketch_params);
        if sketch_params != query_params {
            panic!(
                "Query sketch parameters were not equal to reference sketch parameters. Exiting."
            );
        }
    } else {
        let all_query_range_set =
            (0..command_params.query_files.len()).collect::<FxHashSet<usize>>();
        query_sketches =
            file_io::fastx_to_sketches(&command_params.query_files, &sketch_params, true);
    }

    info!("Generating sketch time: {}", now.elapsed().as_secs_f32());
    let now = Instant::now();
    if command_params.mode == params::Mode::Dist {
        let all_screened_refs;
        let js = (0..query_sketches.len()).collect::<Vec<usize>>();
        let anis: Mutex<Vec<_>> = Mutex::new(vec![]);
        let index_to_screened_refs: Mutex<Vec<_>> =
            Mutex::new(vec![FxHashSet::default(); query_sketches.len()]);
        if !command_params.screen {
            all_screened_refs = (0..ref_sketches.len()).collect::<FxHashSet<usize>>();
        }
        //Get refs to actually sketch
        else {
            info!("Screening to find references to seed");
            let temp_screened_refs: Mutex<FxHashSet<usize>> = Mutex::new(FxHashSet::default());
            js.into_par_iter().for_each(|j| {
                let query_sketch = &query_sketches[j];
                let screened_refs = screen::screen_refs(
                    command_params.screen_val,
                    &kmer_to_sketch,
                    query_sketch,
                    &sketch_params,
                    &ref_sketches,
                );
                debug!(
                    "{} has {} refs passing screening.",
                    query_sketch.file_name,
                    screened_refs.len()
                );
                let mut locked = temp_screened_refs.lock().unwrap();
                locked.extend(&screened_refs);
                let mut locked = index_to_screened_refs.lock().unwrap();
                locked[j] = screened_refs;
            });
            all_screened_refs = temp_screened_refs.into_inner().unwrap();
        }
        let index_to_screened_refs = index_to_screened_refs.into_inner().unwrap();
        //Sketch refs
        let ref_sketches_past_screen;
        if command_params.screen {
            info!(
                "{} references passed screening. Seeding ...",
                all_screened_refs.len()
            );

            let mut temp = ref_sketches.clone();
            file_io::seed_screened_sketches(&sketch_params, &mut temp, &all_screened_refs);
            ref_sketches_past_screen = temp;
        } else {
            ref_sketches_past_screen = ref_sketches;
        }

        info!("Seeding complete");
        //Now compute everything
        let js = (0..query_sketches.len()).collect::<Vec<usize>>();
        js.into_par_iter().for_each(|j| {
            let query_sketch = &query_sketches[j];
            let screened_refs;
            if command_params.screen {
                screened_refs = &index_to_screened_refs[j];
            } else {
                screened_refs = &all_screened_refs;
            }
            screened_refs.into_par_iter().for_each(|i| {
                let ref_sketch = &ref_sketches_past_screen[*i];
                let map_params = chain::map_params_from_sketch(
                    ref_sketch,
                    sketch_params.use_aa,
                    &command_params,
                );
                let ani_res;
                if map_params != params::MapParams::default() {
                    ani_res = chain::chain_seeds(ref_sketch, &query_sketch, map_params);
                } else {
                    ani_res = types::AniEstResult::default();
                }
                let mut locked = anis.lock().unwrap();
                locked.push(ani_res);
            });
        });
        let anis = anis.into_inner().unwrap();
        file_io::write_query_ref_list(&anis, &command_params.out_file_name);
    } else if command_params.mode == params::Mode::Triangle {
        let anis: Mutex<FxHashMap<usize, FxHashMap<usize, types::AniEstResult>>> =
            Mutex::new(FxHashMap::default());

        if ref_sketches.len() == 0 {
            panic!("No reference fastas/sketches found.")
        }

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
                        &ref_sketches,
                    );
                    debug!(
                        "{} has {} refs passing screening.",
                        ref_sketch_i.file_name,
                        screened_refs.len()
                    );
                } else {
                    screened_refs = (i + 1..ref_sketches.len()).collect::<FxHashSet<usize>>();
                }

                screened_refs.into_par_iter().for_each(|j| {
                    if j > i {
                        let map_params = chain::map_params_from_sketch(
                            ref_sketch_i,
                            sketch_params.use_aa,
                            &command_params,
                        );
                        let ref_sketch_j = &ref_sketches[j];
                        let ani_res = chain::chain_seeds(ref_sketch_i, ref_sketch_j, map_params);
                        let mut locked = anis.lock().unwrap();
                        let mapi = locked.entry(i).or_insert(FxHashMap::default());
                        mapi.insert(j, ani_res);
                    }
                });
            });
        let anis = anis.into_inner().unwrap();
        if command_params.sparse {
            file_io::write_sparse_matrix(&anis, &ref_sketches, &command_params.out_file_name);
        } else {
            file_io::write_phyllip_matrix(&anis, &ref_sketches, &command_params.out_file_name);
        }
    }
    info!("Alignment time: {}", now.elapsed().as_secs_f32());
}
