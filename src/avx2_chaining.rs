//! AVX2 implementation of the banded anchor-chaining DP.
//!
//! The chaining recurrence in `chain::chain_anchors_ani` is
//!
//! ```text
//! score[i] = max(0, max over j in W(i) of (s(i,j) + score[j]))
//! ```
//!
//! where `W(i)` is a band of predecessors. The band is
//! `BP_CHAIN_BAND / c` anchors wide - 20 at the default `c = 125` - so each
//! `i` has only ~20 candidate predecessors, all of them independent of one
//! another. That is a natural fit for integer SIMD: one vector evaluates 8
//! candidates at a time and the best is picked with a max-reduction.
//!
//! Two properties of the scalar loop make the vector form exact rather than
//! approximate:
//!
//! 1. Anchors within a chunk share a query contig and are sorted by query
//!    position, so the scalar `break` condition
//!    `query_pos[i] - query_pos[j] > past_chain_length || i - j > band`
//!    is monotone in `j`. The visited set is therefore exactly the index
//!    window `[lo(i), i-1]`, with `lo(i)` non-decreasing in `i` - a
//!    two-pointer walk. (The `ref_contig` mismatch `continue` sits before the
//!    `break` in the scalar loop, but because the window test is monotone it
//!    cannot extend the range, only punch holes in it, which here becomes a
//!    per-lane mask.)
//! 2. The scalar loop scans `j` downwards and updates on a strict `>`, so the
//!    winner is the maximum score, ties broken towards the *largest* `j`. The
//!    reduction below applies that same lexicographic `(score, j)` rule, so it
//!    selects the same predecessor and not merely an equally-scoring one.
//!
//! Scores are `i32` here. That is not a precision trade: `score_anchors` is
//! integer-valued by construction, so this is bit-identical to the f64 form.

use crate::params::*;
use std::arch::x86_64::*;

/// Evaluate the chaining DP for one chunk, writing `score_vec` and
/// `pointer_vec`.
///
/// Inputs are the chunk's anchors in struct-of-arrays form: `qpos`, `rpos`,
/// `rctg` and `rev` (0 or 1) must all have the same length, and every position
/// must fit in `i32` (the caller checks this). `score_vec` and `pointer_vec`
/// must be at least that long.
///
/// # Safety
///
/// Requires AVX2. The slice lengths above must hold.
#[target_feature(enable = "avx2")]
pub unsafe fn chain_chunk_avx2(
    qpos: &[i32],
    rpos: &[i32],
    rctg: &[i32],
    rev: &[i32],
    past_chain_length: i32,
    band: usize,
    anchor_score: i32,
    max_gap_length: i32,
    score_vec: &mut [i32],
    pointer_vec: &mut [usize],
) {
    let n = qpos.len();
    debug_assert!(rpos.len() == n && rctg.len() == n && rev.len() == n);
    debug_assert!(score_vec.len() >= n && pointer_vec.len() >= n);

    let zero = _mm256_setzero_si256();
    let ones = _mm256_cmpeq_epi32(zero, zero);
    let v_lin = _mm256_set1_epi32(D_MAX_LIN_LENGTH);
    let v_maxgap = _mm256_set1_epi32(max_gap_length);
    let v_anchor = _mm256_set1_epi32(anchor_score);
    let v_sentinel = _mm256_set1_epi32(i32::MIN);
    let lane_id = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

    // Lower edge of the band; monotone in i, so this only ever moves forward.
    let mut lo: usize = 0;

    for i in 0..n {
        if i > band {
            let bound = i - band;
            if lo < bound {
                lo = bound;
            }
        }
        while lo < i && qpos[i] - qpos[lo] > past_chain_length {
            lo += 1;
        }

        let mut best_score: i32 = 0;
        let mut best_prev_index: usize = i;

        if lo < i {
            let v_qi = _mm256_set1_epi32(qpos[i]);
            let v_ri = _mm256_set1_epi32(rpos[i]);
            let v_ci = _mm256_set1_epi32(rctg[i]);
            let v_vi = _mm256_set1_epi32(rev[i]);
            let reverse_i = rev[i] != 0;

            // Running best over the band, held per lane and folded at the end.
            let mut acc_score = v_sentinel;
            let mut acc_j = _mm256_set1_epi32(-1);

            let mut j = lo;
            while j < i {
                let rem = i - j;
                let in_range = _mm256_cmpgt_epi32(_mm256_set1_epi32(rem as i32), lane_id);
                let v_j = _mm256_add_epi32(_mm256_set1_epi32(j as i32), lane_id);

                // Full 8-wide loads for interior blocks; the tail block uses
                // masked loads so it never reads past the chunk.
                let (q_j, r_j, c_j, rv_j, s_j) = if rem >= 8 {
                    (
                        _mm256_loadu_si256(qpos.as_ptr().add(j) as *const __m256i),
                        _mm256_loadu_si256(rpos.as_ptr().add(j) as *const __m256i),
                        _mm256_loadu_si256(rctg.as_ptr().add(j) as *const __m256i),
                        _mm256_loadu_si256(rev.as_ptr().add(j) as *const __m256i),
                        _mm256_loadu_si256(score_vec.as_ptr().add(j) as *const __m256i),
                    )
                } else {
                    (
                        _mm256_maskload_epi32(qpos.as_ptr().add(j), in_range),
                        _mm256_maskload_epi32(rpos.as_ptr().add(j), in_range),
                        _mm256_maskload_epi32(rctg.as_ptr().add(j), in_range),
                        _mm256_maskload_epi32(rev.as_ptr().add(j), in_range),
                        _mm256_maskload_epi32(score_vec.as_ptr().add(j), in_range),
                    )
                };

                let d_q = _mm256_sub_epi32(v_qi, q_j);
                let d_r = if reverse_i {
                    _mm256_sub_epi32(r_j, v_ri)
                } else {
                    _mm256_sub_epi32(v_ri, r_j)
                };

                // Everything score_anchors rejects, as lane masks.
                let m_ctg = _mm256_cmpeq_epi32(c_j, v_ci);
                let m_rev = _mm256_cmpeq_epi32(rv_j, v_vi);
                let m_rne = _mm256_xor_si256(_mm256_cmpeq_epi32(r_j, v_ri), ones);
                let m_qne = _mm256_xor_si256(_mm256_cmpeq_epi32(q_j, v_qi), ones);
                let m_dq = _mm256_xor_si256(_mm256_cmpgt_epi32(d_q, v_lin), ones);
                let m_dr = _mm256_xor_si256(_mm256_cmpgt_epi32(d_r, v_lin), ones);
                let m_drpos = _mm256_cmpgt_epi32(d_r, zero);

                let pre = _mm256_and_si256(
                    _mm256_and_si256(m_ctg, m_rev),
                    _mm256_and_si256(m_rne, m_qne),
                );
                let pre = _mm256_and_si256(
                    pre,
                    _mm256_and_si256(_mm256_and_si256(m_dq, m_dr), m_drpos),
                );

                // Zero the offsets on already-rejected lanes before the gap
                // subtraction. Those lanes are masked out below regardless;
                // this just keeps `d_r - d_q` from overflowing i32 on a lane
                // whose offsets are arbitrarily large.
                let d_q_s = _mm256_and_si256(d_q, pre);
                let d_r_s = _mm256_and_si256(d_r, pre);
                let gap = _mm256_abs_epi32(_mm256_sub_epi32(d_r_s, d_q_s));
                let m_gap = _mm256_xor_si256(_mm256_cmpgt_epi32(gap, v_maxgap), ones);

                let valid = _mm256_and_si256(_mm256_and_si256(pre, m_gap), in_range);
                let cand = _mm256_add_epi32(_mm256_sub_epi32(v_anchor, gap), s_j);
                let cand = _mm256_blendv_epi8(v_sentinel, cand, valid);

                // (score, j) lexicographic max: strictly better score wins,
                // an equal score keeps the larger j.
                let gt = _mm256_cmpgt_epi32(cand, acc_score);
                let eq = _mm256_cmpeq_epi32(cand, acc_score);
                let j_max = _mm256_max_epi32(acc_j, v_j);
                acc_j = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(acc_j, j_max, eq),
                    v_j,
                    gt,
                );
                acc_score = _mm256_max_epi32(acc_score, cand);

                j += 8;
            }

            let mut lane_scores = [0i32; 8];
            let mut lane_js = [0i32; 8];
            _mm256_storeu_si256(lane_scores.as_mut_ptr() as *mut __m256i, acc_score);
            _mm256_storeu_si256(lane_js.as_mut_ptr() as *mut __m256i, acc_j);

            // Fold the lanes under the same rule. The seed (0, i) reproduces
            // the scalar loop's initial state, so a candidate scoring 0 loses
            // to it exactly as the scalar strict `>` would have it.
            for l in 0..8 {
                if lane_js[l] < 0 {
                    continue;
                }
                let s = lane_scores[l];
                let jj = lane_js[l] as usize;
                if s > best_score || (s == best_score && jj > best_prev_index) {
                    best_score = s;
                    best_prev_index = jj;
                }
            }
        }

        score_vec[i] = best_score;
        pointer_vec[i] = best_prev_index;
    }
}
