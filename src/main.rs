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
use skani::search;
use skani::screen;
use skani::types;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::sync::Mutex;
use std::time::Instant;
use std::path::Path;
fn main() {
    let matches = Command::new("skani")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1")
        .about("skani")
        .subcommand(
            SubCommand::with_name("help").setting(AppSettings::Hidden)
        )
        .subcommand(
            SubCommand::with_name(params::SKETCH_STRING)
            .about("Sketch (index) genomes. Usage: skani sketch genome1.fa genome2.fa ... -o new_sketch_folder")
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
                        .help("file with each line containing one fasta/sketch file.")
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
                .arg(Arg::new("output").short('o').help("output folder. Creates a folder if it does not exist, and overwrites the contents in folder if it does.").takes_value(true).required(true))
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
            .about("Compute ANI/AAI for queries against references fasta files or pre-computed sketch files. Usage: skani dist query.fa ref1.fa ref2.fa ... or use -q/--ql and -r/--rl options.")
                .arg(Arg::new("v").short('v').help("verbose output."))
                .arg(Arg::new("trace").long("trace").help("trace level output."))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
                .arg(
                    Arg::new("query")
                        .index(1)
                        .help("query fasta or sketch.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("reference")
                        .index(2)
                        .help("reference fasta(s) or sketch(es).")
                        .takes_value(true)
                        .multiple(true)
                )
                .arg(
                    Arg::new("queries")
                        .short('q')
                        .help("query fasta(s) or sketch(es)")
                        .takes_value(true)
                        .multiple(true)
                )
                .arg(
                    Arg::new("references")
                        .short('r')
                        .help("reference fasta(s) or sketch(es)")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("reference list file")
                        .long("rl")
                        .help("file with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("query list file")
                        .long("ql")
                        .help("file with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .help_heading("OUTPUT")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .help("output file name; rewrites file by default (default: output to stdout)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("n")
                        .short('n')
                        .help("max number of results to show for each query. (default: unlimited)")
                        .takes_value(true)
                )
                .help_heading("ALGORITHM PARAMETERS")
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size (default: ANI-15, AAI-6).")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help("compression factor (default: ANI-80, AAI-15).")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .help("calculate AAI instead (default: ANI)."),
                )
                .group(
                    ArgGroup::new("ref")
                        .arg("reference")
                        .arg("references")
                        .arg("reference list file")
//                        .required(true)
                )
                .group(
                    ArgGroup::new("q")
                        .arg("query")
                        .arg("queries")
                        .arg("query list file")
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
            .about("Compute a lower triangular distance ANI/AAI matrix in .phyllip format. Usage: skani triangle genome1.fa genome2.fa genome3.fa ...")
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
                        .help("file with each line containing one fasta/sketch file.")
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
                .arg(Arg::new("output").short('o').help("output. ").takes_value(true))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .arg(Arg::new("s").short('s').takes_value(true).help("screen. ").hidden(true))
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
        .subcommand(
            SubCommand::with_name(params::SEARCH_STRING)
            .about("Search or screen queries against a pre-sketched database of reference genomes in a memory efficient manner. Usage: skani search -d sketch_folder query1.fa query2.fa ... ")
                .arg(Arg::new("v").short('v').help("verbose output."))
                .arg(Arg::new("trace").long("trace").help("trace level output."))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
                .arg(
                    Arg::new("sketched database folder")
                        .short('d')
                        .help("folder of outputs from `skani sketch`.")
                        .takes_value(true)
                        .required(true)
                )
                .arg(
                    Arg::new("query")
                        .index(1)
                        .help("query fasta(s) or sketch(es).")
                        .multiple(true)
                        .takes_value(true),
                )
                .arg(
                    Arg::new("queries")
                        .short('q')
                        .help("query fasta(s) or sketch(es)")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("query list file")
                        .long("ql")
                        .help("file with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .help_heading("OUTPUT")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .help("output file name; rewrites file by default (default: output to stdout).")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("n")
                        .short('n')
                        .help("max number of results to show for each query. (default: unlimited)")
                        .takes_value(true)
                )
                .group(
                    ArgGroup::new("q")
                        .arg("query")
                        .arg("queries")
                        .arg("query list file")
                        .required(true),
                )
                .help_heading("ALGORITHM PARAMETERS")
                .arg(Arg::new("s").short('s').takes_value(true).help("ANI/AAI screening threshold."))
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
        .get_matches();

    let (mut sketch_params, command_params) = parse::parse_params(&matches);
    if command_params.ref_files.len() == 0 {
        panic!("No reference fastas/sketches found.")
    }

    //SKETCHING
    if command_params.mode == params::Mode::Sketch {
        let p = format!("{}", command_params.out_file_name);
        std::fs::create_dir_all(p).unwrap();

        let num_iters = command_params.ref_files.len();
        (0..num_iters).into_par_iter().for_each(|i| {
            let ref_sketches = file_io::fastx_to_sketches(
                &vec![command_params.ref_files[i].clone()],
                &sketch_params,
                !command_params.screen,
            );
            let marker_ref_sketches = ref_sketches.iter().map(|x| types::Sketch::get_markers_only(x)).collect::<Vec<types::Sketch>>();
            if ref_sketches.len() > 0{
                let sketch = &ref_sketches[0];
                let marker_sketch = &marker_ref_sketches[0];
                let path = Path::new(&sketch.file_name);
                let filename = path.file_name().unwrap().to_str().unwrap();
                let mut file_bin = BufWriter::new(
                    File::create(format!("{}/{}.sketch", &command_params.out_file_name, filename)).unwrap(),
                );
                let mut file_bin_marker = BufWriter::new(
                    File::create(format!("{}/{}.marker", &command_params.out_file_name, filename)).unwrap(),
                );

                bincode::serialize_into(
                    &mut file_bin,
                    &(
                        &sketch_params,
                        sketch,
                    ),
                )
                .unwrap();

                bincode::serialize_into(
                    &mut file_bin_marker,
                    &(
                        &sketch_params,
                        marker_sketch,
                    ),
                )
                .unwrap();
            }

        });
        return;
    }

    if command_params.mode == params::Mode::Search{
        search::search(&command_params);
        return;
    }
    //

    let now = Instant::now();
    let ref_sketches;
    let kmer_to_sketch;

    if command_params.refs_are_sketch {
        info!("Sketches detected.");
        (sketch_params, ref_sketches) =
            file_io::sketches_from_sketch(&command_params.ref_files, command_params.mode == params::Mode::Search);
        if command_params.screen {
            kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches);
        } else {
            kmer_to_sketch = types::KmerToSketch::default();
        }
    } else {
        if command_params.screen && command_params.mode == params::Mode::Dist {
            ref_sketches =
                file_io::fastx_to_sketches(&command_params.ref_files, &sketch_params, false);
        } else {
            ref_sketches =
                file_io::fastx_to_sketches(&command_params.ref_files, &sketch_params, true);
        }
        if command_params.screen {
            kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches)
        } else {
            kmer_to_sketch = types::KmerToSketch::default();
        }
    }

    let query_sketches;
    let query_params;
    if command_params.queries_are_sketch {
        (query_params, query_sketches) =
            file_io::sketches_from_sketch(&command_params.query_files, false);
        if sketch_params != query_params {
            panic!(
                "Query sketch parameters were not equal to reference sketch parameters. Exiting."
            );
        }
    } else {
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
                if ani_res.ani > 0.5 {
                    let mut locked = anis.lock().unwrap();
                    locked.push(ani_res);
                }
            });
        });
        let anis = anis.into_inner().unwrap();
        file_io::write_query_ref_list(
            &anis,
            &command_params.out_file_name,
            command_params.max_results,
        );
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
