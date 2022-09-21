use crate::types::*;
use fxhash::FxHashMap;
pub const D_FRAGMENT_LENGTH: usize = 200000;
pub const D_MAX_GAP_LENGTH: f64 = 50.;
pub const D_ANCHOR_SCORE: f64 = 50.;
pub const D_MIN_ANCHORS: usize = 5;
pub const D_LENGTH_CUTOFF: usize = D_FRAGMENT_LENGTH;
pub const D_FRAC_COVER_CUTOFF: f64 = 0.02;
pub const D_CHAIN_BAND: usize = 100;

#[derive(Default)]
pub struct MapParams {
    pub fragment_length: usize,
    pub max_gap_length: f64,
    pub anchor_score: f64,
    pub min_anchors: usize,
    pub length_cutoff: usize,
    pub mode: String,
    pub frac_cover_cutoff: f64,
    pub chain_band: usize,
    pub k: usize,
}

pub fn fragment_length_formula(n: usize) -> usize {
    return (n as f64).sqrt() as usize * 500;
}

#[derive(Default)]
pub struct SketchParams {
    pub cs: Vec<usize>,
    pub ks: Vec<usize>,
    pub use_syncs: bool,
    pub use_aa: bool,
    pub acgt_to_aa_encoding: Vec<KmerBits>,
}

impl SketchParams {
    pub fn new(cs: Vec<usize>, ks: Vec<usize>, use_syncs: bool, use_aa: bool) -> SketchParams {
        let mut acgt_to_aa_encoding = vec![0;64];
        const DNA_TO_AA: [u8; 64] =
            *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
        let letter_to_int_aa: FxHashMap<u8, KmerBits> = [
            (b'A', 0),
            (b'R', 1),
            (b'D', 2),
            (b'C', 3),
            (b'E', 4),
            (b'F', 5),
            (b'G', 6),
            (b'H', 7),
            (b'I', 8),
            (b'K', 9),
            (b'L', 10),
            (b'M', 11),
            (b'P', 12),
            (b'Q', 13),
            (b'R', 14),
            (b'S', 15),
            (b'T', 16),
            (b'V', 17),
            (b'W', 18),
            (b'Y', 19),
            (b'*', 20),
        ]
        .iter()
        .cloned()
        .collect();
        for i in 0..64{
            acgt_to_aa_encoding[i] = letter_to_int_aa[&DNA_TO_AA[i]];
        }
        return SketchParams {
            cs,
            ks,
            use_syncs,
            use_aa,
            acgt_to_aa_encoding,
        };
    }
}
