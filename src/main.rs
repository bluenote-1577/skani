use clap::{AppSettings, Arg, Command, SubCommand};
use skani::types::*;
use skani::chain;
use skani::file_io::fasta_to_sketch;

fn main() {
    let matches = Command::new("skani")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1")
        .about("skani")
        .arg(
                Arg::new("query")
                    .index(1)
                    .help("query fasta.")
                    .takes_value(true)
                    .required(true)
            )
            .arg(
                Arg::new("reference")
                    .index(2)
                    .help("reference fasta.")
                    .takes_value(true)
                    .required(true)
            ).
            arg(
                Arg::new("k")
                    .short('k')
                    .help("k-mer size.")
                    .default_value("30")
                    .takes_value(true)
            ).
            arg(
                Arg::new("c")
                    .short('c')
                    .help("compression factor.")
                    .default_value("100")
                    .takes_value(true)
            ).get_matches();

    let ref_file = matches.value_of("reference").unwrap();
    let query_file = matches.value_of("query").unwrap();
    let k = matches.value_of("k").unwrap().parse::<usize>().unwrap();
    let c = matches.value_of("c").unwrap().parse::<usize>().unwrap();
    let ref_sketch = fasta_to_sketch(ref_file, k, c);
    let query_sketch = fasta_to_sketch(query_file, k, c);

    chain::chain_seeds(&ref_sketch, &query_sketch, k, c);
}
