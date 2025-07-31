use clap::Parser;
use std::env;
use skani::cli::{Cli, Commands};
use skani::dist;
use skani::parse;
use skani::search;
use skani::sketch;
use skani::triangle;

//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though. 
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let cli = Cli::parse();

    let (sketch_params, command_params) = parse::parse_params_from_cli(&cli);

    let cmd_txt = env::args().into_iter().collect::<Vec<String>>().join(" ");
    let log_str = &cmd_txt[0..usize::min(cmd_txt.len(), 250)];
    if cmd_txt.len() > 250{
        log::info!("{} ...", log_str);
    }
    else{
        log::info!("{}", log_str);
    }

    match cli.command {
        Commands::Sketch(_) => {
            sketch::sketch(command_params, sketch_params);
        },
        Commands::Search(_) => {
            search::search(command_params);
        },
        Commands::Dist(_) => {
            dist::dist(command_params, sketch_params);
        },
        Commands::Triangle(_) => {
            triangle::triangle(command_params, sketch_params);
        },
    }
}
