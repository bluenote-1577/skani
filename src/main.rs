use clap::{AppSettings, Arg, ArgGroup, Command, SubCommand};
use std::env;
use skani::dist;
use skani::cmd_line::*;
use skani::params;
use skani::parse;
use skani::search;
use skani::sketch;
use skani::triangle;

//Use this allocator instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Not an issue when not using MUSL though. 
//#[global_allocator]
//static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let matches = Command::new("skani")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.0.1")
        .about("fast, robust ANI calculation and database searching for metagenomic contigs and assemblies. \n\nQuick ANI calculation:\nskani dist genome1.fa genome2.fa \n\nMemory-efficient database search:\nskani sketch genomes/* -o database; skani search -d database query1.fa query2.fa ...\n\nAll-to-all comparison:\nskani triangle genomes/*")
        .subcommand(
            SubCommand::with_name("help").setting(AppSettings::Hidden)
        )
        .subcommand(
            SubCommand::with_name(params::SKETCH_STRING)
            .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("3")
                        .help("Number of threads. ")
                        .takes_value(true),
                )
            .arg(
                    Arg::new("avx2")
                        .long("avx2")
                )
            .about("Sketch (index) genomes. Usage: skani sketch genome1.fa genome2.fa ... -o new_sketch_folder")
                .help_heading("INPUT/OUTPUT")
                .arg(
                    Arg::new("fasta_files")
                        .index(1)
                        .help("fastas to sketch.")
                        .takes_value(true)
                        .multiple(true)
                )
                .arg(
                    Arg::new("fasta_list")
                        .short('l')
                        .help("File with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .arg(Arg::new("output sketch folder").short('o').help("Output folder where sketch files are placed. Creates a folder if it does not exist, and overwrites the contents in folder if it does.").takes_value(true).required(true).display_order(1))
                .help_heading("SKETCH PARAMETERS")
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .hidden(true)
                        .help("Use amino acid to calculate AAI instead.\t[default: ANI]"),
                )
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size.\t[default: 15]")
                        .takes_value(true)
                        .hidden(true)
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help(H_C_FACTOR)
                        .takes_value(true),
                )
                .arg(
                    Arg::new("marker_c")
                        .short('m')
                        .help("Marker k-mer compression factor. Markers are used for filtering. You want at least ~100 markers, so genome_size/marker_c > 100 is highly recommended. Higher value is more time/memory efficient. \t[default: 1000] ")
                        .takes_value(true),
                )
                .group(
                    ArgGroup::new("ref")
                        .arg("fasta_files")
                        .arg("fasta_list")
                        .required(true),
                )
                
                .help_heading("MISC")
                .arg(Arg::new("v").short('v').help("Debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("Trace level verbosity."))

            )
        .subcommand(
            SubCommand::with_name(params::DIST_STRING)
            .about("Compute ANI for queries against references fasta files or pre-computed sketch files. Usage: skani dist query.fa ref1.fa ref2.fa ... or use -q/--ql and -r/--rl options.")
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("3")
                        .help("Number of threads.")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .hidden(true)
                        .help("Use amino acid to calculate AAI instead.\t[default: ANI]"),
                )
                .arg(
                    Arg::new("query")
                        .index(1)
                        .help("Query fasta or sketch.")
                        .takes_value(true)
                )
                .arg(
                    Arg::new("reference")
                        .index(2)
                        .help("Reference fasta(s) or sketch(es).")
                        .takes_value(true)
                        .multiple(true)
                )
                .arg(
                    Arg::new("queries")
                        .short('q')
                        .help("Query fasta(s) or sketch(es)")
                        .takes_value(true)
                        .multiple(true)
                )
                .arg(
                    Arg::new("references")
                        .short('r')
                        .help("Reference fasta(s) or sketch(es)")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("reference list file")
                        .long("rl")
                        .help("File with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("query list file")
                        .long("ql")
                        .help("File with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new(IND_CTG_QRY)
                        .long(CMD_IND_CTG_QRY)
                        .help(H_IND_CTG_QRY)
                )
                .arg(
                    Arg::new(IND_CTG_REF)
                        .long(CMD_IND_CTG_REF)
                        .help(H_IND_CTG_REF)
                )
                .help_heading("OUTPUT")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .help("Output file name; rewrites file by default\t[default: output to stdout]")
                        .takes_value(true)
                        .display_order(1)
                )
                .arg(
                    Arg::new(MIN_ALIGN_FRAC)
                        .long(CMD_MIN_ALIGN_FRAC)
                        .help(H_MIN_ALIGN_FRAC)
                        .takes_value(true)
                        .display_order(100)
                )
                .arg(
                    Arg::new("n")
                        .short('n')
                        .help("Max number of results to show for each query.\t[default: unlimited]")
                        .takes_value(true)
                )
                .arg(
                    Arg::new(CONF_INTERVAL)
                        .long(CMD_CONF_INTERVAL)
                        .help(H_CONF_INTERVAL)
                )
                .help_heading("ALGORITHM PARAMETERS")
                .arg(
                    Arg::new("marker_c")
                        .short('m')
                        .help("Marker k-mer compression factor. Markers are used for filtering. You want at least ~100 markers, so genome_size/marker_c > 100 is highly recommended. Higher value is more time/memory efficient. \t[default: 1000] ")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size.\t[default: 15].")
                        .takes_value(true)
                        .hidden(true)
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help(H_C_FACTOR)
                        .takes_value(true),
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
                .arg(Arg::new("s").short('s').takes_value(true).help(H_SCREEN))
                .arg(
                    Arg::new(ROBUST)
                        .long(CMD_ROBUST)
                        .help(H_ROBUST),
                )
                .arg(
                    Arg::new("median")
                        .long("median")
                        .help("Estimate median identity instead of average (mean) identity."),
                )
                .arg(
                    Arg::new(FULL_INDEX)
                        .long(CMD_FULL_INDEX)
                        .help(H_FULL_INDEX),
                )
                .help_heading("MISC")
                .arg(Arg::new("v").short('v').help("Debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("Trace level verbosity."))
        )
        .subcommand(
            SubCommand::with_name(params::TRIANGLE_STRING)
            .about("Compute a lower triangular distance ANI/AF matrix. Usage: skani triangle genome1.fa genome2.fa genome3.fa ...")
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("3")
                        .help("Number of threads.")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
                .arg(
                    Arg::new("fasta_list")
                        .short('l')
                        .help("File with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .help("Use amino acid to calculate AAI instead.\t[default: ANI]"),
                )
                .arg(
                    Arg::new("fasta_files")
                        .index(1)
                        .help("Fasta(s) or sketch(es).")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("individual contig")
                        .short('i')
                        .help("Use individual sequences instead the entire file for multi-fastas.")
                )
                .help_heading("OUTPUT")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .help("Output file name; rewrites file by default\t[default: output to stdout]")
                        .takes_value(true)
                        .display_order(1)
                )
                .arg(
                    Arg::new(FULL_MAT)
                        .long(CMD_FULL_MAT)
                        .help(H_FULL_MAT)
                )
                .arg(
                    Arg::new(MIN_ALIGN_FRAC)
                        .long(CMD_MIN_ALIGN_FRAC)
                        .help(H_MIN_ALIGN_FRAC)
                        .takes_value(true)
                )
                .arg(
                    Arg::new(CONF_INTERVAL)
                        .long(CMD_CONF_INTERVAL)
                        .help(H_CONF_INTERVAL_TRI)
                )
                .arg(
                    Arg::new("sparse")
                        .long("sparse")
                        .short('E')
                        .help("Output comparisons in a row-by-row form (i.e. sparse matrix) in the same form as `skani dist`. Only pairs with aligned fraction > --min-af are output."),
                )
                .help_heading("ALGORITHM PARAMETERS")
                .arg(
                    Arg::new("marker_c")
                        .short('m')
                        .help("Marker k-mer compression factor. Markers are used for filtering. You want at least ~100 markers, so genome_size/marker_c > 100 is highly recommended. Higher value is more time/memory efficient. \t[default: 1000] ")
                        .takes_value(true),
                )
                .arg(Arg::new("s").short('s').takes_value(true).help(H_SCREEN))
                .arg(
                    Arg::new("k")
                        .short('k')
                        .help("k-mer size.\t[default: 15]")
                        .takes_value(true)
                        .hidden(true)
                )
                .arg(
                    Arg::new("c")
                        .short('c')
                        .help(H_C_FACTOR)
                        .takes_value(true),
                )
                .group(
                    ArgGroup::new("ref")
                        .arg("fasta_files")
                        .arg("fasta_list")
                        .required(true),
                )
                .arg(
                    Arg::new(ROBUST)
                        .long(CMD_ROBUST)
                        .help(H_ROBUST),
                )
                .arg(
                    Arg::new("median")
                        .long("median")
                        .help("Estimate median identity instead of average (mean) identity."),
                )
                .help_heading("MISC")
                .arg(Arg::new("v").short('v').help("Debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("Trace level verbosity."))
        )
        .subcommand(
            SubCommand::with_name(params::SEARCH_STRING)
            .about("Search queries against a large pre-sketched database of reference genomes in a memory efficient manner. Algorithm parameters are determined by skani sketch parameters. Usage: skani search -d sketch_folder query1.fa query2.fa ... ")
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("3")
                        .help("Number of threads.")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
                .arg(
                    Arg::new("sketched database folder")
                        .short('d')
                        .help("Output folder from `skani sketch`.")
                        .takes_value(true)
                        .required(true)
                )
                .arg(
                    Arg::new("query")
                        .index(1)
                        .help("Query fasta(s) or sketch(es).")
                        .multiple(true)
                        .takes_value(true),
                )
                .arg(
                    Arg::new("queries")
                        .short('q')
                        .help("Query fasta(s) or sketch(es)")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("query list file")
                        .long("ql")
                        .help("File with each line containing one fasta/sketch file.")
                        .takes_value(true),
                )
                .arg(
                    Arg::new(IND_CTG_QRY)
                        .long(CMD_IND_CTG_QRY)
                        .help(H_IND_CTG_QRY)
                )
                .help_heading("OUTPUT")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .help("Output file name; rewrites file by default\t[default: output to stdout].")
                        .takes_value(true)
                        .display_order(1)
                )
                .arg(
                    Arg::new(CONF_INTERVAL)
                        .long(CMD_CONF_INTERVAL)
                        .help(H_CONF_INTERVAL)
                )
                .arg(
                    Arg::new(MIN_ALIGN_FRAC)
                        .long(CMD_MIN_ALIGN_FRAC)
                        .help(H_MIN_ALIGN_FRAC)
                        .takes_value(true)
                )
                .arg(
                    Arg::new("n")
                        .short('n')
                        .help("Max number of results to show for each query.\t[default: unlimited]")
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
                .arg(
                    Arg::new(KEEP_REFS)
                        .long(CMD_KEEP_REFS)
                        .help(H_KEEP_REFS),
                )
                .arg(
                    Arg::new(FULL_INDEX)
                        .long(CMD_FULL_INDEX)
                        .help(H_FULL_INDEX),
                )
                .arg(Arg::new("s").short('s').takes_value(true).help(H_SCREEN))
                .arg(
                    Arg::new(ROBUST)
                        .long(CMD_ROBUST)
                        .help(H_ROBUST),
                )
                .arg(
                    Arg::new("median")
                        .long("median")
                        .help("Estimate median identity instead of average (mean) identity."),
                )
                .help_heading("MISC")
                .arg(Arg::new("v").short('v').help("Debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("Trace level verbosity."))

        )
        .get_matches();

    let (sketch_params, command_params) = parse::parse_params(&matches);

    let cmd_txt = env::args().into_iter().collect::<Vec<String>>().join(" ");
    log::info!("{}", cmd_txt);

    //SKETCHING
    if command_params.mode == params::Mode::Sketch {
        sketch::sketch(command_params, sketch_params);
    } else if command_params.mode == params::Mode::Search {
        search::search(command_params);
    } else if command_params.mode == params::Mode::Dist {
        dist::dist(command_params, sketch_params);
    } else if command_params.mode == params::Mode::Triangle {
        triangle::triangle(command_params, sketch_params);
    }
}
