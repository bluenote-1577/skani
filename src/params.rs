use crate::types::*;
use gbdt::gradient_boost::GBDT;

pub const VERSION: &str = "0.3.0";

pub const GB_IN_BYTES: usize = 1_073_741_824;
pub const SMALL_VEC_SIZE: usize = 1;
pub const KMER_SK_SMALL_VEC_SIZE: usize = 3;
pub const INTERMEDIATE_WRITE_COUNT: usize = 5000;
//pub const INTERMEDIATE_WRITE_COUNT: usize = 2;
use serde::{Serialize, Deserialize};
use fxhash::FxHashMap;
pub const D_FRAGMENT_LENGTH: usize = 200000;
pub const STOP_CODON: MarkerBits = 21;
pub const DEFAULT_C: &str = "125";
pub const DEFAULT_C_AAI: &str = "15";
pub const DEFAULT_K: &str = "15";
pub const DEFAULT_K_AAI: &str = "6";
pub const D_MAX_GAP_LENGTH: f64 = 300.;
pub const D_MAX_GAP_LENGTH_AAI: f64 = 50.;
pub const D_MAX_LIN_LENGTH: f64 = 5000.;
pub const D_ANCHOR_SCORE_ANI: f64 = 20.;
pub const D_ANCHOR_SCORE_AAI: f64 = 20.;
pub const D_MIN_ANCHORS_ANI: usize = 3;
pub const D_MIN_ANCHORS_AAI: usize = 5;
pub const D_LENGTH_CUTOFF: usize = D_FRAGMENT_LENGTH;
pub const D_FRAC_COVER_CUTOFF: &str = "15";
pub const D_ANI_AND_COVER_CUTOFF: f64 = 0.95;
pub const D_FRAC_COVER_CUTOFF_AA: &str = "5";
//pub const D_CHAIN_BAND: usize = 50;
//pub const D_CHAIN_BAND_AAI: usize = 125;
pub const ORF_SIZE: usize = 30;
pub const MARKER_C_DEFAULT: &str = "1000";
pub const K_MARKER_AA: usize = 10;
pub const K_MARKER_DNA: usize = 21;
pub const SEARCH_STRING: &str = "search";
pub const DIST_STRING: &str = "dist";
pub const SKETCH_STRING: &str = "sketch";
pub const TRIANGLE_STRING: &str = "triangle";
pub const CHUNK_SIZE_DNA: usize = 20000;
pub const CHUNK_SIZE_AA: usize = 20000;
pub const MIN_LENGTH_CONTIG: usize = 500;
pub const MIN_LENGTH_COVER_AAI: usize = 500;
pub const MIN_LENGTH_COVER: usize = 500;
pub const BP_CHAIN_BAND: usize = 2500;
pub const BP_CHAIN_BAND_AAI: usize = 500;
pub const SEARCH_AAI_CUTOFF_DEFAULT: f64 = 0.60;
pub const SEARCH_ANI_CUTOFF_DEFAULT: f64 = 0.80;
pub const SCREEN_MINIMUM_KMERS: usize = 20;
pub const FULL_INDEX_THRESH: usize = 50;
pub const REPET_KMER_THRESHOLD: usize = 8_000_000;
pub const OVERLAP_ORTHOLOGOUS_FRACTION: f32  = 0.50;
pub const TOTAL_BASES_REGRESS_CUTOFF: usize = 150000;
pub const LEARNED_INFO_HELP: &str = "Learned ANI mode detected. ANI may be adjusted according to a regression model trained on MAGs.";

pub const FAST_C: usize = 200;
pub const SLOW_C: usize = 30;
pub const MEDIUM_C: usize = 70;
pub const SMALL_M: usize = 200;

pub const ASCII_N: usize = 78;
pub const ASCII_N_SMALL: usize = 110;



#[derive(PartialEq)]
pub enum Mode {
    Sketch,
    Dist,
    Triangle,
    Search,
}

#[derive(Default)]
pub struct MapParams<'a> {
    pub fragment_length: usize,
    pub max_gap_length: f64,
    pub anchor_score: f64,
    pub min_anchors: usize,
    pub length_cutoff: usize,
    pub frac_cover_cutoff: f64,
    pub length_cover_cutoff: usize,
    pub index_chain_band: usize,
    pub k: usize,
    pub amino_acid: bool,
    pub min_score: f64,
    pub robust: bool,
    pub median: bool,
    pub bp_chain_band: usize,
    pub min_length_cover: usize,
    pub model: Option<&'a GBDT>
}

#[derive(PartialEq)]
pub struct CommandParams{
    pub screen: bool,
    pub screen_val: f64,
    pub mode: Mode,
    pub out_file_name: String,
    pub ref_files: Vec<String>,
    pub query_files: Vec<String>,
    pub refs_are_sketch: bool,
    pub queries_are_sketch: bool,
    pub robust: bool,
    pub median: bool,
    pub sparse: bool,
    pub full_matrix: bool,
    pub diagonal: bool,
    pub max_results: usize,
    pub individual_contig_q: bool,
    pub individual_contig_r: bool,
    pub min_aligned_frac: f64,
    pub keep_refs: bool,
    pub est_ci: bool,
    pub learned_ani: bool,
    pub detailed_out: bool,
    pub distance: bool,
    pub rescue_small: bool,
    pub separate_sketches: bool,
}

pub fn fragment_length_formula(_n: usize, aa: bool) -> usize {
//    return (n as f64).sqrt() as usize * 10;
    if aa{
        CHUNK_SIZE_AA
    }
    else{
        CHUNK_SIZE_DNA
    }
//    return (n as f64).sqrt() as usize * 3;
}

#[derive(Default,  PartialEq, Serialize, Deserialize, Debug, Clone)]
pub struct SketchParams {
    pub c: usize,
    pub k: usize,
    pub marker_c: usize,
    pub use_syncs: bool,
    pub use_aa: bool,
    pub acgt_to_aa_encoding: Vec<MarkerBits>,
    pub acgt_to_aa_letters: Vec<u8>,
    pub orf_size: usize,
}

impl SketchParams {
    pub fn new(marker_c: usize, c: usize, k: usize, use_syncs: bool, use_aa: bool) -> SketchParams {
        let mut acgt_to_aa_encoding = vec![0;64];
                let letter_to_int_aa: FxHashMap<u8, MarkerBits> = [
            (b'A', 0),
            (b'R', 1),
            (b'N', 2),
            (b'D', 3),
            (b'C', 4),
            (b'E', 5),
            (b'F', 6),
            (b'G', 7),
            (b'H', 8),
            (b'I', 9),
            (b'K', 10),
            (b'L', 11),
            (b'M', 12),
            (b'P', 13),
            (b'Q', 14),
            (b'R', 15),
            (b'S', 16),
            (b'T', 17),
            (b'V', 18),
            (b'W', 19),
            (b'Y', 20),
            (b'*', STOP_CODON),
        ]
        .iter()
        .cloned()
        .collect();
        for i in 0..64{
            acgt_to_aa_encoding[i] = letter_to_int_aa[&DNA_TO_AA[i]];
        }
        let orf_size = ORF_SIZE;
        let marker_c = marker_c;
        if c > marker_c{
            panic!("We currently don't allow c ({}) > m ({}). -m should be larger than c.", c,  marker_c);
        }
        SketchParams {
            c,
            k,
            marker_c,
            use_syncs,
            use_aa,
            acgt_to_aa_encoding,
            acgt_to_aa_letters: DNA_TO_AA.to_vec(),
            orf_size,
        }
    }
}
