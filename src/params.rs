pub const D_FRAGMENT_LENGTH: usize = 200000;
pub const D_MAX_GAP_LENGTH: f64 = 50.;
pub const D_ANCHOR_SCORE: f64 = 50.;
pub const D_MIN_ANCHORS: usize = 5;
pub const D_LENGTH_CUTOFF: usize = D_FRAGMENT_LENGTH;
pub const D_FRAC_COVER_CUTOFF: f64 = 0.02;
pub const D_CHAIN_BAND: usize = 100;

#[derive(Default)]
pub struct MapParams{
    pub fragment_length: usize,
    pub max_gap_length: f64,
    pub anchor_score: f64,
    pub min_anchors: usize,
    pub length_cutoff : usize,
    pub mode : String,
    pub frac_cover_cutoff: f64,
    pub chain_band: usize,
    pub k: usize,
}

pub fn fragment_length_formula(n : usize) -> usize{
    return (n as f64).sqrt() as usize * 500;
}

//#[derive(Default)]
//pub struct SketchParams{
//    pub c: usize,
//    pub k: usize,
//    pub guess_params: bool
//}
