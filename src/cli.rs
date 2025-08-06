use clap::{Parser, Subcommand, Args};

#[derive(Parser)]
#[clap(
    name = "skani",
    version = "0.3.0",
    about = "fast, robust ANI calculation and database searching for metagenomic contigs and assemblies. \n\nQuick ANI calculation:\nskani dist genome1.fa genome2.fa \n\nMemory-efficient database search:\nskani sketch genomes/* -o database; skani search -d database query1.fa query2.fa ...\n\nAll-to-all comparison:\nskani triangle genomes/*",
    arg_required_else_help = true, disable_help_subcommand = true
)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Sketch (index) genomes.
    /// Usage: skani sketch genome1.fa genome2.fa ... -o new_sketch_folder
    Sketch(SketchArgs),
    
    /// Compute ANI for queries against references fasta files or pre-computed sketch files.
    /// Usage: skani dist query.fa ref1.fa ref2.fa ... or use -q/--ql and -r/--rl options.
    Dist(DistArgs),
    
    /// Compute a lower triangular ANI/AF matrix.
    /// Usage: skani triangle genome1.fa genome2.fa genome3.fa ...
    Triangle(TriangleArgs),
    
    /// Search queries against a large pre-sketched database of reference genomes in a memory efficient manner.
    /// Usage: skani search -d sketch_folder query1.fa query2.fa ...
    Search(SearchArgs),
}

#[derive(Args)]
#[clap(group(
    clap::ArgGroup::new("input_group")
        .required(true)
))]
pub struct SketchArgs {
    /// Number of threads
    #[clap(short = 't', default_value = "3")]
    pub threads: String,

    /// fastas to sketch
    #[clap(help_heading = "INPUT/OUTPUT", group = "input_group")]
    pub fasta_files: Vec<String>,
    
    /// File with each line containing one fasta/sketch file
    #[clap(short = 'l', help_heading = "INPUT/OUTPUT", group = "input_group")]
    pub fasta_list: Option<String>,
    
    /// Use individual sequences instead the entire file for multi-fastas. CURRENTLY DOES NOT WORK WITH `skani search`.
    #[clap(short = 'i', help_heading = "INPUT/OUTPUT")]
    pub individual_contig: bool,
    
    /// Output folder where sketch files are placed. Creates a folder if it does not exist, and overwrites the contents in folder if it does.
    #[clap(short = 'o', required = true, display_order = 1, help_heading = "INPUT/OUTPUT")]
    pub output: String,

    /// Slower skani mode; 4x slower and more memory. Gives much more accurate AF for distant genomes. More accurate ANI for VERY fragmented assemblies (< 3kb N50), but less accurate ANI otherwise. Alias for -c 30.
    #[clap(long = "slow", help_heading = "PRESETS")]
    pub slow: bool,
    
    /// Medium skani mode; 2x slower and more memory. More accurate AF and more accurate ANI for moderately fragmented assemblies (< 10kb N50). Alias for -c 70.
    #[clap(long = "medium", help_heading = "PRESETS")]
    pub medium: bool,
    
    /// Faster skani mode; 2x faster and less memory. Less accurate AF and less accurate ANI for distant genomes, but works ok for high N50 and > 95% ANI. Alias for -c 200.
    #[clap(long = "fast", help_heading = "PRESETS")]
    pub fast: bool,

    /// Create separate .sketch files instead of consolidated database format
    #[clap(long = "separate-sketches", help_heading = "INPUT/OUTPUT")]
    pub separate_sketches: bool,

    /// Use amino acid to calculate AAI instead. [default: ANI]
    #[clap(short = 'a', long = "aai", hide = true, help_heading = "SKETCH PARAMETERS")]
    pub aai: bool,
    
    /// k-mer size. [default: 15]
    #[clap(short = 'k', hide = true, help_heading = "SKETCH PARAMETERS")]
    pub k: Option<String>,
    
    /// Compression factor (k-mer subsampling rate). [default: 125]
    #[clap(short = 'c', help_heading = "SKETCH PARAMETERS")]
    pub c: Option<String>,
    
    /// Marker k-mer compression factor. Markers are used for filtering. Consider decreasing to ~200-300 if working with small genomes (e.g. plasmids or viruses). [default: 1000]
    #[clap(short = 'm', help_heading = "SKETCH PARAMETERS")]
    pub marker_c: Option<String>,

    /// Debug level verbosity
    #[clap(short = 'v', long = "debug", help_heading = "MISC")]
    pub debug: bool,
    
    /// Trace level verbosity
    #[clap(long = "trace", help_heading = "MISC")]
    pub trace: bool,
}

#[derive(Args)]
#[clap(group(
    clap::ArgGroup::new("query_group")
        .required(true)
))]
pub struct DistArgs {
    /// Number of threads
    #[clap(short = 't', default_value = "3")]
    pub threads: String,

    /// Use amino acid to calculate AAI instead. [default: ANI]
    #[clap(short = 'a', long = "aai", hide = true, help_heading = "INPUTS")]
    pub aai: bool,
    
    /// Query fasta or sketch
    #[clap(help_heading = "INPUTS", group = "query_group")]
    pub query: Option<String>,
    
    /// Reference fasta(s) or sketch(es)
    #[clap(help_heading = "INPUTS")]
    pub reference: Vec<String>,
    
    /// Query fasta(s) or sketch(es)
    #[clap(short = 'q', multiple_values = true, help_heading = "INPUTS", group = "query_group")]
    pub queries: Vec<String>,
    
    /// Reference fasta(s) or sketch(es)
    #[clap(short = 'r', multiple_values = true, help_heading = "INPUTS")]
    pub references: Vec<String>,
    
    /// File with each line containing one fasta/sketch file
    #[clap(long = "rl", help_heading = "INPUTS")]
    pub reference_list: Option<String>,
    
    /// File with each line containing one fasta/sketch file
    #[clap(long = "ql", help_heading = "INPUTS", group = "query_group")]
    pub query_list: Option<String>,
    
    /// Use individual sequences for the QUERY in a multi-line fasta
    #[clap(long = "qi", help_heading = "INPUTS")]
    pub qi: bool,
    
    /// Use individual sequences for the REFERENCE in a multi-line fasta
    #[clap(long = "ri", help_heading = "INPUTS")]
    pub ri: bool,

    /// Output file name; rewrites file by default [default: output to stdout]
    #[clap(short = 'o', display_order = 1, help_heading = "OUTPUT")]
    pub output: Option<String>,
    
    /// Only output ANI values where one genome has aligned fraction > than this value. [default: 15]
    #[clap(long = "min-af", display_order = 100, help_heading = "OUTPUT")]
    pub min_af: Option<String>,
    
    /// Max number of results to show for each query. [default: unlimited]
    #[clap(short = 'n', help_heading = "OUTPUT")]
    pub n: Option<String>,
    
    /// Output [5%,95%] ANI confidence intervals using percentile bootstrap on the putative ANI distribution
    #[clap(long = "ci", help_heading = "OUTPUT")]
    pub ci: bool,
    
    /// Print additional info including contig N50s and more
    #[clap(long = "detailed", help_heading = "OUTPUT")]
    pub detailed: bool,
    
    /// Only display the first part of contig names (before first whitespace)
    #[clap(long = "short-header", help_heading = "OUTPUT")]
    pub short_header: bool,

    /// Slower skani mode; 4x slower and more memory. Gives much more accurate AF for distant genomes. More accurate ANI for VERY fragmented assemblies (< 3kb N50), but less accurate ANI otherwise. Alias for -c 30.
    #[clap(long = "slow", help_heading = "PRESETS")]
    pub slow: bool,
    
    /// Medium skani mode; 2x slower and more memory. More accurate AF and more accurate ANI for moderately fragmented assemblies (< 10kb N50). Alias for -c 70.
    #[clap(long = "medium", help_heading = "PRESETS")]
    pub medium: bool,
    
    /// Faster skani mode; 2x faster and less memory. Less accurate AF and less accurate ANI for distant genomes, but works ok for high N50 and > 95% ANI. Alias for -c 200.
    #[clap(long = "fast", help_heading = "PRESETS")]
    pub fast: bool,
    
    /// Mode for small genomes such as viruses or plasmids (< 20 kb). Can be much faster for large data, but is slower/less accurate on bacterial-sized genomes. Alias for: -c 30 -m 200 --faster-small.
    #[clap(long = "small-genomes", help_heading = "PRESETS")]
    pub small_genomes: bool,

    /// Disable regression model for ANI prediction. [default: learned ANI used for c >= 70 and >= 150,000 bases aligned and not on individual contigs]
    #[clap(long = "no-learned-ani", help_heading = "ALGORITHM PARAMETERS")]
    pub no_learned_ani: bool,
    
    /// Marker k-mer compression factor. Markers are used for filtering. Consider decreasing to ~200-300 if working with small genomes (e.g. plasmids or viruses). [default: 1000]
    #[clap(short = 'm', help_heading = "ALGORITHM PARAMETERS")]
    pub marker_c: Option<String>,
    
    /// k-mer size. [default: 15]
    #[clap(short = 'k', hide = true, help_heading = "ALGORITHM PARAMETERS")]
    pub k: Option<String>,
    
    /// Compression factor (k-mer subsampling rate). [default: 125]
    #[clap(short = 'c', help_heading = "ALGORITHM PARAMETERS")]
    pub c: Option<String>,
    
    /// Screen out pairs with *approximately* < % identity using k-mer sketching. [default: 80]
    #[clap(short = 's', help_heading = "ALGORITHM PARAMETERS")]
    pub s: Option<String>,
    
    /// Estimate mean after trimming off 10%/90% quantiles
    #[clap(long = "robust", help_heading = "ALGORITHM PARAMETERS")]
    pub robust: bool,
    
    /// Estimate median identity instead of average (mean) identity
    #[clap(long = "median", help_heading = "ALGORITHM PARAMETERS")]
    pub median: bool,
    
    /// Do not use hash-table inverted index for faster ANI filtering. [default: load index if > 100 query files or using the --qi option]
    #[clap(long = "no-marker-index", help_heading = "ALGORITHM PARAMETERS")]
    pub no_marker_index: bool,
    
    /// Filter genomes with < 20 marker k-mers more aggressively. Much faster for many small genomes but may miss some comparisons.
    #[clap(long = "faster-small", help_heading = "ALGORITHM PARAMETERS")]
    pub faster_small: bool,

    /// Debug level verbosity
    #[clap(short = 'v', long = "debug", help_heading = "MISC")]
    pub debug: bool,
    
    /// Trace level verbosity
    #[clap(long = "trace", help_heading = "MISC")]
    pub trace: bool,
}

#[derive(Args)]
#[clap(group(
    clap::ArgGroup::new("input_group")
        .required(true)
))]
pub struct TriangleArgs {
    /// Number of threads
    #[clap(short = 't', default_value = "3")]
    pub threads: String,

    /// File with each line containing one fasta/sketch file
    #[clap(short = 'l', help_heading = "INPUTS", group = "input_group")]
    pub fasta_list: Option<String>,
    
    /// Use amino acid to calculate AAI instead. [default: ANI]
    #[clap(short = 'a', long = "aai", hide = true, help_heading = "INPUTS")]
    pub aai: bool,
    
    /// Fasta(s) or sketch(es)
    #[clap(help_heading = "INPUTS", group = "input_group")]
    pub fasta_files: Vec<String>,
    
    /// Use individual sequences instead the entire file for multi-fastas
    #[clap(short = 'i', help_heading = "INPUTS")]
    pub individual_contig: bool,

    /// Output file name; rewrites file by default [default: output to stdout]
    #[clap(short = 'o', display_order = 1, help_heading = "OUTPUT")]
    pub output: Option<String>,
    
    /// Output full matrix instead of lower-triangular matrix
    #[clap(long = "full-matrix", help_heading = "OUTPUT")]
    pub full_matrix: bool,
    
    /// Output the diagonal of the ANI matrix (i.e. self-self comparisons) for both dense and sparse matrices
    #[clap(long = "diagonal", help_heading = "OUTPUT")]
    pub diagonal: bool,
    
    /// Only output ANI values where one genome has aligned fraction > than this value. [default: 15]
    #[clap(long = "min-af", help_heading = "OUTPUT")]
    pub min_af: Option<String>,
    
    /// Output [5%,95%] ANI confidence intervals using percentile bootstrap on the putative ANI distribution. Only works with --sparse or -E.
    #[clap(long = "ci", help_heading = "OUTPUT")]
    pub ci: bool,
    
    /// Print additional info including contig N50s and more
    #[clap(long = "detailed", help_heading = "OUTPUT")]
    pub detailed: bool,
    
    /// Only display the first part of contig names (before first whitespace)
    #[clap(long = "short-header", help_heading = "OUTPUT")]
    pub short_header: bool,
    
    /// Output 100 - ANI instead of ANI, creating a distance instead of a similarity matrix. No effect if using --sparse or -E.
    #[clap(long = "distance", help_heading = "OUTPUT")]
    pub distance: bool,
    
    /// Output comparisons in a row-by-row form (i.e. sparse matrix) in the same form as `skani dist`. Only pairs with aligned fraction > --min-af are output.
    #[clap(long = "sparse", short = 'E', help_heading = "OUTPUT")]
    pub sparse: bool,

    /// Slower skani mode; 4x slower and more memory. Gives much more accurate AF for distant genomes. More accurate ANI for VERY fragmented assemblies (< 3kb N50), but less accurate ANI otherwise. Alias for -c 30.
    #[clap(long = "slow", help_heading = "PRESETS")]
    pub slow: bool,
    
    /// Medium skani mode; 2x slower and more memory. More accurate AF and more accurate ANI for moderately fragmented assemblies (< 10kb N50). Alias for -c 70.
    #[clap(long = "medium", help_heading = "PRESETS")]
    pub medium: bool,
    
    /// Faster skani mode; 2x faster and less memory. Less accurate AF and less accurate ANI for distant genomes, but works ok for high N50 and > 95% ANI. Alias for -c 200.
    #[clap(long = "fast", help_heading = "PRESETS")]
    pub fast: bool,
    
    /// Mode for small genomes such as viruses or plasmids (< 20 kb). Can be much faster for large data, but is slower/less accurate on bacterial-sized genomes. Alias for: -c 30 -m 200 --faster-small.
    #[clap(long = "small-genomes", help_heading = "PRESETS")]
    pub small_genomes: bool,

    /// Disable regression model for ANI prediction. [default: learned ANI used for c >= 70 and >= 150,000 bases aligned and not on individual contigs]
    #[clap(long = "no-learned-ani", help_heading = "ALGORITHM PARAMETERS")]
    pub no_learned_ani: bool,
    
    /// Marker k-mer compression factor. Markers are used for filtering. Consider decreasing to ~200-300 if working with small genomes (e.g. plasmids or viruses). [default: 1000]
    #[clap(short = 'm', help_heading = "ALGORITHM PARAMETERS")]
    pub marker_c: Option<String>,
    
    /// Screen out pairs with *approximately* < % identity using k-mer sketching. [default: 80]
    #[clap(short = 's', help_heading = "ALGORITHM PARAMETERS")]
    pub s: Option<String>,
    
    /// k-mer size. [default: 15]
    #[clap(short = 'k', hide = true, help_heading = "ALGORITHM PARAMETERS")]
    pub k: Option<String>,
    
    /// Compression factor (k-mer subsampling rate). [default: 125]
    #[clap(short = 'c', help_heading = "ALGORITHM PARAMETERS")]
    pub c: Option<String>,
    
    /// Estimate mean after trimming off 10%/90% quantiles
    #[clap(long = "robust", help_heading = "ALGORITHM PARAMETERS")]
    pub robust: bool,
    
    /// Estimate median identity instead of average (mean) identity
    #[clap(long = "median", help_heading = "ALGORITHM PARAMETERS")]
    pub median: bool,
    
    /// Filter genomes with < 20 marker k-mers more aggressively. Much faster for many small genomes but may miss some comparisons.
    #[clap(long = "faster-small", help_heading = "ALGORITHM PARAMETERS")]
    pub faster_small: bool,

    /// Debug level verbosity
    #[clap(short = 'v', long = "debug", help_heading = "MISC")]
    pub debug: bool,
    
    /// Trace level verbosity
    #[clap(long = "trace", help_heading = "MISC")]
    pub trace: bool,
}

#[derive(Args)]
#[clap(group(
    clap::ArgGroup::new("query_group")
        .required(true)
))]
pub struct SearchArgs {
    /// Number of threads
    #[clap(short = 't', default_value = "3")]
    pub threads: String,

    /// Output folder from `skani sketch`
    #[clap(short = 'd', required = true, help_heading = "INPUTS")]
    pub database: String,
    
    /// Query fasta(s) or sketch(es)
    #[clap(multiple_values = true, help_heading = "INPUTS", group = "query_group")]
    pub query: Vec<String>,
    
    /// Query fasta(s) or sketch(es)
    #[clap(short = 'q', multiple_values = true, help_heading = "INPUTS", group = "query_group")]
    pub queries: Vec<String>,
    
    /// File with each line containing one fasta/sketch file
    #[clap(long = "ql", help_heading = "INPUTS", group = "query_group")]
    pub query_list: Option<String>,
    
    /// Use individual sequences for the QUERY in a multi-line fasta
    #[clap(long = "qi", help_heading = "INPUTS")]
    pub qi: bool,

    /// Output file name; rewrites file by default [default: output to stdout]
    #[clap(short = 'o', display_order = 1, help_heading = "OUTPUT")]
    pub output: Option<String>,
    
    /// Output [5%,95%] ANI confidence intervals using percentile bootstrap on the putative ANI distribution
    #[clap(long = "ci", help_heading = "OUTPUT")]
    pub ci: bool,
    
    /// Print additional info including contig N50s and more
    #[clap(long = "detailed", help_heading = "OUTPUT")]
    pub detailed: bool,
    
    /// Only display the first part of contig names (before first whitespace)
    #[clap(long = "short-header", help_heading = "OUTPUT")]
    pub short_header: bool,
    
    /// Only output ANI values where one genome has aligned fraction > than this value. [default: 15]
    #[clap(long = "min-af", help_heading = "OUTPUT")]
    pub min_af: Option<String>,
    
    /// Max number of results to show for each query. [default: unlimited]
    #[clap(short = 'n', help_heading = "OUTPUT")]
    pub n: Option<String>,

    /// Disable regression model for ANI prediction. [default: learned ANI used for c >= 70 and >= 150,000 bases aligned and not on individual contigs]
    #[clap(long = "no-learned-ani", help_heading = "ALGORITHM PARAMETERS")]
    pub no_learned_ani: bool,
    
    /// Keep reference sketches in memory if the sketch passes the marker filter. Takes more memory but is much faster when querying many similar sequences.
    #[clap(long = "keep-refs", help_heading = "ALGORITHM PARAMETERS")]
    pub keep_refs: bool,
    
    /// Do not use hash-table inverted index for faster ANI filtering. [default: load index if > 100 query files or using the --qi option]
    #[clap(long = "no-marker-index", help_heading = "ALGORITHM PARAMETERS")]
    pub no_marker_index: bool,
    
    /// Screen out pairs with *approximately* < % identity using k-mer sketching. [default: 80]
    #[clap(short = 's', help_heading = "ALGORITHM PARAMETERS")]
    pub s: Option<String>,
    
    /// Estimate mean after trimming off 10%/90% quantiles
    #[clap(long = "robust", help_heading = "ALGORITHM PARAMETERS")]
    pub robust: bool,
    
    /// Estimate median identity instead of average (mean) identity
    #[clap(long = "median", help_heading = "ALGORITHM PARAMETERS")]
    pub median: bool,

    /// Debug level verbosity
    #[clap(short = 'v', long = "debug", help_heading = "MISC")]
    pub debug: bool,
    
    /// Trace level verbosity
    #[clap(long = "trace", help_heading = "MISC")]
    pub trace: bool,
}
