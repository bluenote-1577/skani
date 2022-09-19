pub const FRAGMENT_LENGTH: usize = 200000;
pub const MAX_GAP_LENGTH: f64 = 50.;
pub const ANCHOR_SCORE: f64 = 50.;
pub const MIN_ANCHORS: usize = 5;
pub const TRIANGLE_LENGTH_CUTOFF: usize = FRAGMENT_LENGTH;
pub const FRAC_COVER_CUTOFF: f64 = 0.02;

pub struct MapParams{
    pub fragment_length: usize,
    pub max_gap_length: usize,
    pub anchor_score: usize,
    pub min_anchors: usize,
    pub triangle_length_cutoff: usize,
    pub frac_cover_cutoff: f64
}

pub struct SketchParams{
    pub c: usize,
    pub k: usize,
    pub guess_params: bool
}
