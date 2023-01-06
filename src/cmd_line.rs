pub const MIN_ALIGN_FRAC: &str = "min aligned frac";
pub const CMD_MIN_ALIGN_FRAC: &str = "min-aligned-fraction";
pub const H_MIN_ALIGN_FRAC: &str = "Only output ANI/AAi values with aligned fraction > than this value.\t[default: ANI 0.15, AAI 0.05]";

pub const IND_CTG_QRY: &str = "individual contig query";
pub const CMD_IND_CTG_QRY: &str = "qi";
pub const H_IND_CTG_QRY: &str = "Use individual sequences for the QUERY in a multi-line fasta.";

pub const IND_CTG_REF: &str = "individual contig ref";
pub const CMD_IND_CTG_REF: &str = "ri";
pub const H_IND_CTG_REF: &str = "Use individual sequences for the REFERENCE in a multi-line fasta.";

pub const FULL_INDEX: &str = "marker index";
pub const CMD_FULL_INDEX: &str = "marker-index";
pub const H_FULL_INDEX: &str = "Loads a larger index for faster filtering. Uses more memory and takes longer to load, but faster for many all-to-all comparisons.\t[default: load index if > 50 query files, don't load otherwise]";

pub const ROBUST: &str = "robust";
pub const CMD_ROBUST: &str = "robust";
pub const H_ROBUST: &str = "Robust ani estimation; estimate mean after trim off 10%/90% quantiles.";

pub const FULL_MAT: &str = "full-matrix";
pub const CMD_FULL_MAT: &str = "full-matrix";
pub const H_FULL_MAT: &str = "Output full matrix instead of lower-triangular matrix.";

pub const KEEP_REFS: &str = "keep-refs";
pub const CMD_KEEP_REFS: &str = "keep-refs";
pub const H_KEEP_REFS: &str = "Keep reference sketches in memory if the sketch passes the marker filter. Takes more memory but is much faster when querying many similar sequences.";

pub const C_FACTOR: &str = "c";
pub const CMD_C_FACTOR: &str = "c";
pub const H_C_FACTOR: &str = "Compression factor. Memory usage and runtime is inversely proportional to c. Lower c allows for ANI comparison of more distant genomes.\t[default: 120]";

pub const H_SCREEN: &str = "Screen out pairs with < % identity using a hash table in constant time.\t[default: 0.80]";
