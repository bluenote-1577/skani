use crate::params::*;
use crate::cmd_line::*;
use clap::parser::ArgMatches;
use log::LevelFilter;
use log::*;
use std::fs;
use std::fs::File;
use std::io::{prelude::*, BufReader};

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
        Some(SEARCH_STRING) => {
            mode = Mode::Search;
            matches_subc = matches.subcommand_matches(SEARCH_STRING).unwrap();
            //            return parse_params_search(matches_subc);
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

    simple_logging::log_to_stderr(LevelFilter::Info);
    if matches_subc.is_present("v") {
        simple_logging::log_to_stderr(LevelFilter::Debug);
    }
    if matches_subc.is_present("trace") {
        simple_logging::log_to_stderr(LevelFilter::Trace);
    }

    if mode == Mode::Search {
        return parse_params_search(matches_subc);
    }

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
    } else if mode == Mode::Dist {
        if let Some(values) = matches_subc.values_of("reference") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("reference list file") {
            ref_files = vec![];
            ref_file_list = Some(values);
        } else if let Some(values) = matches_subc.values_of("references") {
            ref_files = values.map(|x| x.to_string()).collect();
        } else {
            panic!("No reference inputs found.");
        }
    } else {
        panic!("PATH TODO");
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
    let mut max_results = usize::MAX;
    if mode == Mode::Dist {
        max_results = matches_subc
            .value_of("n")
            .unwrap_or("1000000000")
            .parse::<usize>()
            .unwrap();
        if let Some(values) = matches_subc.values_of("query") {
            query_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.values_of("queries") {
            query_files = values.map(|x| x.to_string()).collect();
        } else if let Some(values) = matches_subc.value_of("query list file") {
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

    let min_aligned_frac;
    if mode != Mode::Sketch{
        let def_maf = if amino_acid{ D_FRAC_COVER_CUTOFF_AA} else {D_FRAC_COVER_CUTOFF};
        min_aligned_frac = matches_subc.value_of(MIN_ALIGN_FRAC).unwrap_or(def_maf).parse::<f64>().unwrap();
    }
    else{
        min_aligned_frac = 0.;
    }

    let marker_c = matches_subc.value_of("marker_c").unwrap_or(MARKER_C).parse::<usize>().unwrap();
    let out_file_name;
    if mode == Mode::Triangle {
        out_file_name = matches_subc.value_of("output").unwrap_or("").to_string();
    } else if mode == Mode::Sketch {
        out_file_name = matches_subc.value_of("output sketch folder").unwrap_or("").to_string();
    } else if mode == Mode::Dist {
        out_file_name = matches_subc.value_of("output").unwrap_or("").to_string();
    }
    else{
        panic!("Mode doesn't exist");
    }

    let mut screen_val = 0.;
    let mut robust = false;
    let mut median = false;
    if mode == Mode::Triangle || mode == Mode::Dist {
        screen_val = matches_subc.value_of("s").unwrap_or("0.00").parse::<f64>().unwrap() / 100.;
    }
    if mode == Mode::Triangle || mode == Mode::Search || mode == Mode::Dist {
        robust = matches_subc.is_present("robust");
        median = matches_subc.is_present("median");
    }

    let sketch_params = SketchParams::new(marker_c, c, k, use_syncs, amino_acid);

    let mut refs_are_sketch = ref_files.len() > 0;
    for ref_file in ref_files.iter() {
        if !ref_file.contains(".sketch") && !ref_file.contains(".marker") && !ref_file.contains("markers.bin"){
            refs_are_sketch = false;
            break;
        }
    }

    let mut queries_are_sketch = query_files.len() > 0;
    for query_file in query_files.iter() {
        if !query_file.contains(".sketch") && !query_file.contains("markers.bin") {
            queries_are_sketch = false;
            break;
        }
    }

    
    let individual_contig_q;
    let individual_contig_r;

    if mode == Mode::Triangle{
        if matches_subc.is_present("individual contig"){
            individual_contig_q = true;
            individual_contig_r = true;
        }
        else{
            individual_contig_q = false;
            individual_contig_r = false;
        }
    }
    else if mode == Mode::Dist{
        individual_contig_q = matches_subc.is_present("individual contig query");
        individual_contig_r = matches_subc.is_present("individual contig ref");
    }
    else{
        individual_contig_q = false;
        individual_contig_r = false;
    }
    
    let full_matrix;
    if mode == Mode::Triangle{
        full_matrix = matches_subc.is_present(FULL_MAT);
    }
    else{
        full_matrix = false;
    }

    let screen;
    if mode == Mode::Dist{
        if query_files.len() > FULL_INDEX_THRESH || individual_contig_q{
            screen = true;
        }
        else{
            screen = matches_subc.is_present(FULL_INDEX);
        }
    }
    else if mode == Mode::Triangle{
        screen = true;
    }
    else{
        screen = false;
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
        full_matrix,
        max_results,
        individual_contig_q,
        individual_contig_r,
        min_aligned_frac,
        keep_refs: false
    };

    return (sketch_params, command_params);
}

pub fn parse_params_search(matches_subc: &ArgMatches) -> (SketchParams, CommandParams) {
    let mode = Mode::Search;
    let out_file_name = matches_subc.value_of("output").unwrap_or("").to_string();

    let mut query_files = vec![];
    let mut query_file_list = None;
    let max_results = matches_subc
        .value_of("n")
        .unwrap_or("1000000000")
        .parse::<usize>()
        .unwrap();
    if let Some(values) = matches_subc.values_of("query") {
        query_files = values.map(|x| x.to_string()).collect();
    } else if let Some(values) = matches_subc.values_of("queries") {
        query_files = values.map(|x| x.to_string()).collect();
    } else if let Some(values) = matches_subc.value_of("query list file") {
        query_file_list = Some(values);
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
    let ref_folder = matches_subc.value_of("sketched database folder").unwrap();
    let paths =
        fs::read_dir(ref_folder).expect("Issue with folder specified by -d option; exiting");
    let ref_files = paths
        .into_iter()
        .map(|x| x.unwrap().path().to_str().unwrap().to_string())
        .collect();
    let refs_are_sketch = true;

    let mut queries_are_sketch = query_files.len() > 0;
    for query_file in query_files.iter() {
        if !query_file.contains(".sketch") && !query_file.contains("markers.bin") {
            queries_are_sketch = false;
            break;
        }
    }

    let robust = matches_subc.is_present("robust");
    let median = matches_subc.is_present("median");
    let sparse = false;

    let screen_val = matches_subc
        .value_of("s")
        .unwrap_or("0.00")
        .parse::<f64>()
        .unwrap();
    let screen;
    let individual_contig_q = matches_subc.is_present("individual contig query");
    if query_files.len() > FULL_INDEX_THRESH || individual_contig_q{
        screen = true;
    }
    else{
        screen = matches_subc.is_present(FULL_INDEX);
    }

    let min_aligned_frac = matches_subc.value_of(MIN_ALIGN_FRAC).unwrap_or("-1.0").parse::<f64>().unwrap();
    let keep_refs = matches_subc.is_present(KEEP_REFS);

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
        full_matrix: false,
        max_results,
        individual_contig_q,
        individual_contig_r: false,
        min_aligned_frac,
        keep_refs,
    };

    if command_params.ref_files.len() == 0 {
        error!("No valid reference fastas or sketches found.");
        std::process::exit(1)
    }

    return (SketchParams::default(), command_params);
}
