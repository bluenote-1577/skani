use crate::params::*;
use clap::parser::ArgMatches;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use log::LevelFilter;
use log::{debug, info};
use rayon::prelude::*;

pub fn parse_params(matches: &ArgMatches) -> (SketchParams, CommandParams) {
    let mode;
    let matches_subc;
    match matches.subcommand_name() {
        Some(SKETCH_STRING) => {
            mode = Mode::Sketch;
            matches_subc = matches.subcommand_matches(SKETCH_STRING).unwrap()
        }
        Some(DIST_STRING) => {
            mode = Mode::Dist;
            matches_subc = matches.subcommand_matches(DIST_STRING).unwrap()
        }
        Some(TRIANGLE_STRING) => {
            mode = Mode::Triangle;
            matches_subc = matches.subcommand_matches(TRIANGLE_STRING).unwrap()
        }
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
    let robust = matches_subc.is_present("robust");
    if matches_subc.is_present("aai") {
        amino_acid = true;
    } else {
        amino_acid = false;
    }

    let mut ref_files: Vec<String>;
    let mut ref_file_list = None;
    if mode == Mode::Triangle {
        if let Some(values) = matches_subc.values_of("fasta_files") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("fasta_list") {
            ref_files = vec![];
            ref_file_list = Some(values);
        } else {
            panic!("No reference inputs found.");
        }
    } else {
        if let Some(values) = matches_subc.values_of("reference") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else {
            panic!("No reference inputs found.");
        }
    }

    if !ref_file_list.is_none(){
        let ref_file_list = ref_file_list.unwrap();
        let file = File::open(ref_file_list).unwrap();
        let reader = BufReader::new(file);
        let mut temp_vec = vec![];

        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        ref_files = temp_vec;
    }

    let mut query_files = vec![];
    if mode == Mode::Dist {
        if let Some(values) = matches_subc.values_of("query") {
            query_files = values.map(|x| x.to_string()).collect();
        } else {
            query_files = vec![];
        }
    }
    let def_k = if amino_acid { DEFAULT_K_AAI } else { DEFAULT_K };
    let def_c = if amino_acid { DEFAULT_C_AAI } else { DEFAULT_C };
    let k = matches_subc
        .value_of("k")
        .unwrap_or(def_k)
        .parse::<usize>()
        .unwrap();
    let ks = vec![k];
    let c;
    let use_syncs = false;
    c = matches_subc
        .value_of("c")
        .unwrap_or(def_c)
        .parse::<usize>()
        .unwrap();
    let cs = vec![c];

    let mut out_file_name = matches_subc
        .value_of("o")
        .unwrap_or("skani_res")
        .to_string();
    if mode == Mode::Triangle {
        out_file_name = format!("{}.matrix", out_file_name).to_string();
    } else if mode == Mode::Sketch {
        out_file_name = format!("{}.sketch", out_file_name).to_string();
    } else if mode == Mode::Dist {
        out_file_name = format!("{}.dist", out_file_name).to_string();
    }
    simple_logging::log_to_stderr(LevelFilter::Info);
    if matches_subc.is_present("v") {
        simple_logging::log_to_stderr(LevelFilter::Debug);
    }
    if matches_subc.is_present("trace") {
        simple_logging::log_to_stderr(LevelFilter::Trace);
    }

    let mut screen = false;
    let mut screen_val = 0.;
    if mode != Mode::Sketch {
        screen = matches_subc.is_present("s");
        if screen {
            screen_val = matches_subc.value_of("s").unwrap().parse::<f64>().unwrap();
        } else {
            screen_val = 0.;
        }
    }

    let sketch_params = SketchParams::new(cs, ks, use_syncs, amino_acid);

    let ref_is_sketch = if ref_files.len() == 1 && ref_files[0].contains(".sketch") {
        true
    } else {
        false
    };
    let query_is_sketch = if query_files.len() == 1 && query_files[0].contains(".sketch") {
        true
    } else {
        false
    };

    let command_params = CommandParams {
        screen,
        screen_val,
        mode,
        out_file_name,
        ref_files,
        query_files,
        ref_is_sketch,
        query_is_sketch,
        robust
    };

    return (sketch_params, command_params);
}
