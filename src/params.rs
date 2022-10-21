use crate::types::*;
pub const SMALL_VEC_SIZE: usize = 1;
use serde::{Serialize, Deserialize};
use fxhash::FxHashMap;
pub const D_FRAGMENT_LENGTH: usize = 200000;
pub const STOP_CODON: KmerBits = 21;
pub const DEFAULT_C: &str = "100";
pub const DEFAULT_C_AAI: &str = "15";
pub const DEFAULT_K: &str = "15";
pub const DEFAULT_K_AAI: &str = "6";
pub const D_MAX_GAP_LENGTH: f64 = 300.;
pub const D_MAX_LIN_LENGTH: f64 = 5000.;
pub const D_ANCHOR_SCORE: f64 = 50.;
pub const D_MIN_ANCHORS: usize = 5;
pub const D_LENGTH_CUTOFF: usize = D_FRAGMENT_LENGTH;
pub const D_FRAC_COVER_CUTOFF: f64 = 0.10;
pub const D_FRAC_COVER_CUTOFF_AA: f64 = 0.02;
pub const D_CHAIN_BAND: usize = 100;
pub const ORF_SIZE: usize = 100;
pub const MARKER_C: usize = 1000;
pub const K_MARKER_AA: usize = 10;
pub const K_MARKER_DNA: usize = 21;
pub const DIST_STRING: &str = "dist";
pub const SKETCH_STRING: &str = "sketch";
pub const TRIANGLE_STRING: &str = "triangle";
pub const CHUNK_SIZE_DNA: usize = 20000;
pub const MIN_LENGTH_CONTIG: usize = 500;

#[derive(PartialEq)]
pub enum Mode {
    Sketch,
    Dist,
    Triangle,
}

#[derive(Default, PartialEq)]
pub struct MapParams {
    pub fragment_length: usize,
    pub max_gap_length: f64,
    pub anchor_score: f64,
    pub min_anchors: usize,
    pub length_cutoff: usize,
    pub frac_cover_cutoff: f64,
    pub length_cover_cutoff: usize,
    pub chain_band: usize,
    pub k: usize,
    pub amino_acid: bool,
    pub min_score: f64,
    pub robust: bool
}

pub struct CommandParams{
    pub screen: bool,
    pub screen_val: f64,
    pub mode: Mode,
    pub out_file_name: String,
    pub ref_files: Vec<String>,
    pub query_files: Vec<String>,
    pub ref_is_sketch: bool,
    pub query_is_sketch: bool,
    pub robust: bool
}

pub fn fragment_length_formula(n: usize, aa: bool) -> usize {
//    return (n as f64).sqrt() as usize * 10;
    if aa{
        return 20000;
    }
    else{
        return CHUNK_SIZE_DNA;
    }
//    return (n as f64).sqrt() as usize * 3;
}

#[derive(Default,  PartialEq, Serialize, Deserialize)]
pub struct SketchParams {
    pub cs: Vec<usize>,
    pub ks: Vec<usize>,
    pub marker_c: usize,
    pub use_syncs: bool,
    pub use_aa: bool,
    pub acgt_to_aa_encoding: Vec<KmerBits>,
    pub acgt_to_aa_letters: Vec<u8>,
    pub orf_size: usize,
}

impl SketchParams {
    pub fn new(cs: Vec<usize>, ks: Vec<usize>, use_syncs: bool, use_aa: bool) -> SketchParams {
        let mut acgt_to_aa_encoding = vec![0;64];
        let DNA_TO_AA: [u8; 64] =
            *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
        let letter_to_int_aa: FxHashMap<u8, KmerBits> = [
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
        let marker_c = MARKER_C;
        return SketchParams {
            cs,
            ks,
            marker_c,
            use_syncs,
            use_aa,
            acgt_to_aa_encoding,
            acgt_to_aa_letters: DNA_TO_AA.to_vec(),
            orf_size
        };
    }
}
