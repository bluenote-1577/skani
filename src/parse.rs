use crate::params::*;
use clap::parser::ArgMatches;
use log::LevelFilter;
use log::{debug, info};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

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
    if matches_subc.is_present("aai") {
        amino_acid = true;
    } else {
        amino_acid = false;
    }

    let mut ref_files: Vec<String>;
    let mut ref_file_list = None;
    let mut sparse = false;
    if mode == Mode::Triangle {
        sparse = matches_subc.is_present("sparse");
    }
    if mode == Mode::Triangle || mode == Mode::Sketch {
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
        } else if let Some(values) = matches_subc.value_of("ref_fasta_list") {
            ref_files = vec![];
            ref_file_list = Some(values);
        } else {
            panic!("No reference inputs found.");
        }
    }

    if !ref_file_list.is_none() {
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
    let mut query_file_list = None;
    if mode == Mode::Dist {
        if let Some(values) = matches_subc.values_of("query") {
            query_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("query_fasta_list") {
            query_file_list = Some(values)
        }
    }

    if !query_file_list.is_none() {
        let query_file_list = query_file_list.unwrap();
        let file = File::open(query_file_list).unwrap();
        let reader = BufReader::new(file);
        let mut temp_vec = vec![];

        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        query_files = temp_vec;
    }

    let def_k = if amino_acid { DEFAULT_K_AAI } else { DEFAULT_K };
    let def_c = if amino_acid { DEFAULT_C_AAI } else { DEFAULT_C };
    let k = matches_subc
        .value_of("k")
        .unwrap_or(def_k)
        .parse::<usize>()
        .unwrap();
    let c;
    let use_syncs = false;
    c = matches_subc
        .value_of("c")
        .unwrap_or(def_c)
        .parse::<usize>()
        .unwrap();

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

    let mut screen_val = 0.;
    let mut robust = false;
    let mut median = false;
    let screen = matches_subc.is_present("s");
    if mode != Mode::Sketch {
        if screen {
            screen_val = matches_subc.value_of("s").unwrap().parse::<f64>().unwrap();
        } else {
            screen_val = 0.;
        }
        robust = matches_subc.is_present("robust");
        median = matches_subc.is_present("median");
    }

    let sketch_params = SketchParams::new(c, k, use_syncs, amino_acid, !screen && mode == Mode::Triangle);

    let mut refs_are_sketch = ref_files.len() > 0;
    for ref_file in ref_files.iter() {
        if !ref_file.contains(".sketch") {
            refs_are_sketch = false;
            break;
        }
    }

    let mut queries_are_sketch = query_files.len() > 0;
    for query_file in query_files.iter() {
        if !query_file.contains(".sketch") {
            queries_are_sketch = false;
            break;
        }
    }

    let command_params = CommandParams {
        screen,
        screen_val,
        mode,
        out_file_name,
        ref_files,
        query_files,
        refs_are_sketch,
        queries_are_sketch,
        robust,
        median,
        sparse,
    };

    return (sketch_params, command_params);
}
