use clap::{AppSettings, Arg, Command, SubCommand};
use rayon::prelude::*;
use std::sync::Mutex;
use log::LevelFilter;
use log::{info, debug};
use simple_logging;
use skani::chain;
use skani::screen;
use skani::file_io;
use skani::params;
use skani::types;
use std::time::Instant;
fn main() {
    let dist = "dist";
    let classify = "classify";
    let matches = Command::new("skani")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1")
        .about("skani")
        .subcommand(
            SubCommand::with_name(classify)
                .arg(
                    Arg::new("query")
                        .short('q')
                        .help("query fasta.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .help("reference fasta.")
                        .takes_value(true)
                        .multiple(true),
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
                    Arg::new("d")
                        .short('d')
                        .help("dowsample factor for open syncmers.")
                        .takes_value(true),
                )
                .arg(Arg::new("v").short('v').help("verbose.")),
        )
        .subcommand(
            SubCommand::with_name("dist")
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
                    Arg::new("d")
                        .short('d')
                        .help("dowsample factor for open syncmers. ")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("query")
                        .short('q')
                        .help("query fasta.")
                        .takes_value(true),
                )
                .arg(Arg::new("aai").short('a').long("aai").help("use amino acid alphabet."))
                .arg(Arg::new("v").short('v').help("verbose."))
                .arg(Arg::new("trace").long("trace").help("verbose."))
                .arg(Arg::new("m").short('m').help("distance matrix. "))
                .arg(Arg::new("e").short('e').help("eukaryotic option. Changes mapping parameters"))
                .arg(Arg::new("o").short('o').help("output. ").takes_value(true))
                .arg(Arg::new("t").short('t').default_value("20").help("threads. ").takes_value(true))
                .arg(Arg::new("s").short('s').takes_value(true).help("screen. "))
        )
        .get_matches();

    let mode;
    let matches_subc;
    match matches.subcommand_name() {
        Some("classify") => {
            mode = classify;
            matches_subc = matches.subcommand_matches(classify).unwrap()
        } // clone was used
        Some("dist") => {
            mode = dist;
            matches_subc = matches.subcommand_matches(dist).unwrap()
        } // push was used
        _ => {
            panic!()
        } // Either no subcommand or one not tested for...
    }
    let threads = matches_subc.value_of("t").unwrap();
    let threads = threads.parse::<usize>().unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let amino_acid;
    if matches_subc.is_present("aai") {
        amino_acid = true;
    } else {
        amino_acid = false;
    }

    let ref_files: Vec<&str>;
    if let Some(values) = matches_subc.values_of("reference") {
        ref_files = values.collect();
    } else {
        ref_files = vec![];
    }
    let def_k = if amino_acid{params::DEFAULT_K_AAI} else {params::DEFAULT_K};
    let def_c = if amino_acid{params::DEFAULT_C_AAI} else {params::DEFAULT_C};
    let k = matches_subc
        .value_of("k")
        .unwrap_or(def_k)
        .parse::<usize>()
        .unwrap();
    let ks = vec![k];
    let c;
    let use_syncs;
    if !matches_subc.value_of("d").is_none() {
        use_syncs = true;
        c = matches_subc.value_of("d").unwrap().parse().unwrap();
    } else {
        use_syncs = false;
        c = matches_subc
            .value_of("c")
            .unwrap_or(def_c)
            .parse::<usize>()
            .unwrap();
    }
    let cs = vec![c];
    
    let euk = matches_subc.is_present("e");
    let mat_file_name = matches_subc.value_of("o").unwrap_or("skani_res");
    simple_logging::log_to_stderr(LevelFilter::Info);
//    if mode == dist && matches_subc.is_present("m") {
//        simple_logging::log_to_stderr(LevelFilter::Warn);
//    } else {
//        simple_logging::log_to_stderr(LevelFilter::Info);
//    }
    if matches_subc.is_present("v") {
        simple_logging::log_to_stderr(LevelFilter::Debug);
    }
    if matches_subc.is_present("trace") {
        simple_logging::log_to_stderr(LevelFilter::Trace);
    }
    let screen = matches_subc.is_present("s");
    let screen_val;
    if screen{
        screen_val = matches_subc.value_of("s").unwrap().parse::<f64>().unwrap();
    }
    else{
        screen_val = 0.;
    }

    let sketch_params = params::SketchParams::new(cs, ks, use_syncs, amino_acid, euk);

    let now = Instant::now();
    let ref_sketches = file_io::fastx_to_sketches(ref_files, &sketch_params);
    let kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches);
    let query_sketches;
    let triangle;
    if mode == classify {
        triangle = false;
        let query_file = matches_subc.value_of("query").unwrap();
        query_sketches = file_io::fastx_to_multiple_sketch_rewrite(query_file, &sketch_params);
    } else if mode == dist {
        let query_file = matches_subc.value_of("query");
        if query_file.is_none() {
            triangle = true;
            query_sketches = vec![];
        } else {
            triangle = false;
            query_sketches = file_io::fastx_to_sketches(
                vec![query_file.unwrap()],
                &sketch_params,
            );
        }
    } else {
        panic!("MODE ERROR");
    }
    info!("Generating sketch time: {}", now.elapsed().as_secs_f32());

    
    let now = Instant::now();
    if mode == classify {
        for query_sketch in query_sketches.iter() {
            let now = Instant::now();
            for ref_sketch in ref_sketches.iter() {
                let map_params = chain::map_params_from_sketch(ref_sketch, mode, amino_acid, euk, screen);
                chain::chain_seeds(ref_sketch, query_sketch, map_params);
            }
            info!("Alignment time: {}", now.elapsed().as_secs_f32());
        }
    } else if mode == dist {
        if triangle {
            let anis :Mutex<Vec<_>> =
                Mutex::new(vec![vec![types::AniEstResult::default(); ref_sketches.len()]; ref_sketches.len()]);
            (0..ref_sketches.len()-1).collect::<Vec<usize>>().into_par_iter().for_each( |i| {
                let screened_refs;
                let ref_sketch_i = &ref_sketches[i];
                if screen{
                    screened_refs = screen::screen_refs(screen_val, &kmer_to_sketch, &ref_sketch_i, &sketch_params);
                    info!("{} has {} refs passing screening.", ref_sketch_i.file_name, screened_refs.len());
                }
                else{
                    screened_refs = (i+1..ref_sketches.len()).collect::<Vec<usize>>();
                }
                for j in i+1..ref_sketches.len(){
//                    {
//                        let mut locked = anis.lock().unwrap();
//                        locked[i][j].ani = 0.;
//                    }
                    if screened_refs.contains(&j){
                        let map_params = chain::map_params_from_sketch(ref_sketch_i, mode, amino_acid, euk, screen);
                        let ref_sketch_j = &ref_sketches[j];
                        let ani_res = chain::chain_seeds(ref_sketch_i, ref_sketch_j, map_params);
                        let mut locked = anis.lock().unwrap();
                        locked[i][j] = ani_res;
                    }
                }
            });
            let anis = anis.into_inner().unwrap();
            file_io::write_phyllip_matrix(&anis, &ref_sketches, &mat_file_name);
        }
        else{
            let anis: Mutex<Vec<_>>  = Mutex::new(vec![]);
            let screened_refs;
            if screen{
                screened_refs = screen::screen_refs(screen_val, &kmer_to_sketch, &query_sketches[0], &sketch_params);
                info!("{} has {} refs passing screening.", &query_sketches[0].file_name, screened_refs.len());
            }
            else{
                screened_refs = (0..ref_sketches.len()).collect::<Vec<usize>>();
            }
            screened_refs.into_par_iter().for_each( |i| {
                let ref_sketch = &ref_sketches[i];
                let map_params = chain::map_params_from_sketch(ref_sketch, mode, amino_acid, euk, screen);
                let ani_res = chain::chain_seeds(ref_sketch, &query_sketches[0], map_params);
                let mut locked = anis.lock().unwrap();
                locked.push(ani_res);
            });
            let anis = anis.into_inner().unwrap();
            file_io::write_query_ref_list(&anis, &mat_file_name);
        }
    }
    info!("Alignment time: {}", now.elapsed().as_secs_f32());
}
