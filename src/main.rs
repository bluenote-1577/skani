use clap::{AppSettings, Arg, ArgGroup, Command, SubCommand};
use skani::dist;
use skani::params;
use skani::parse;
use skani::search;
use skani::sketch;
use skani::triangle;
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
                .arg(Arg::new("v").short('v').help("debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("trace level verbosity.")),
        )
        .subcommand(
            SubCommand::with_name(params::DIST_STRING)
            .about("Compute ANI/AAI for queries against references fasta files or pre-computed sketch files. Usage: skani dist query.fa ref1.fa ref2.fa ... or use -q/--ql and -r/--rl options.")
                .arg(Arg::new("v").short('v').help("debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("trace level verbosity."))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
                .arg(
                    Arg::new("aai")
                        .short('a')
                        .long("aai")
                        .help("calculate AAI instead (default: ANI)."),
                )
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
                .arg(
                    Arg::new("individual contig query")
                        .long("qi")
                        .help("Calculate ANI on individual sequences for the QUERY instead the entire file for multi-fastas.")
                )
                .arg(
                    Arg::new("individual contig ref")
                        .long("ri")
                        .help("Calculate ANI on individual sequences for the REFERENCE instead the entire file for multi-fastas.")
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
                .arg(Arg::new("s").short('s').takes_value(true).help("screen out pairs with < % identity using a hash table in constant time. (default ANI : 75%, AAI: 50%)"))
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
                .arg(Arg::new("v").short('v').help("debug level verbosity."))
                .arg(Arg::new("trace").long("trace").help("trace level verbosity."))
                .arg(
                    Arg::new("t")
                        .short('t')
                        .default_value("20")
                        .help("threads. ")
                        .takes_value(true),
                )
                .help_heading("INPUTS")
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
                .arg(
                    Arg::new("fasta_files")
                        .index(1)
                        .help("fasta(s) or sketch(es).")
                        .takes_value(true)
                        .multiple(true),
                )
                .arg(
                    Arg::new("individual contig")
                        .short('i')
                        .help("Calculate ANI on individual sequences instead the entire file for multi-fastas.")
                )
                .help_heading("OUTPUT")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .help("output file name; rewrites file by default (default: output to stdout)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("sparse")
                        .long("sparse")
                        .help("output sparse matrix for only non-zero ANI/AAI in the same form as `skani dist`."),
                )
                .help_heading("ALGORITHM PARAMETERS")
                .arg(Arg::new("s").short('s').takes_value(true).help("screen out pairs with < % identity using a hash table in constant time. (default ANI : 75%, AAI: 50%)"))
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
                ,
        )
        .subcommand(
            SubCommand::with_name(params::SEARCH_STRING)
            .about("Search queries against a pre-sketched database of reference genomes in a memory efficient manner. Usage: skani search -d sketch_folder query1.fa query2.fa ... ")
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
                .arg(
                    Arg::new("individual contig query")
                        .long("qi")
                        .help("Calculate ANI on individual sequences for the QUERY instead the entire file for multi-fastas.")
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
                .arg(Arg::new("s").short('s').takes_value(true).help("screen out pairs with < % identity. (default ANI : 75%, AAI: 50%)"))
                .arg(
                    Arg::new("robust")
                        .long("robust")
                        .help("Robust identity estimation; trim off 10%/90% quantiles."),
                )
                .arg(
                    Arg::new("median")
                        .long("median")
                        .help("Estimate median identity instead of average (mean)identity."),
                ),
        )
        .get_matches();

    let (sketch_params, command_params) = parse::parse_params(&matches);

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
