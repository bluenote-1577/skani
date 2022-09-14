use crate::types::*;
use std::ops::Bound::*;
use bio::data_structures::interval_tree::IntervalTree;
use partitions::*;

pub fn score_anchors(anchor_curr: &Anchor, anchor_past: &Anchor) -> f64 {
    if anchor_curr.ref_contig != anchor_past.query_contig {
        return f64::MIN;
    }
    if anchor_curr.reverse_match != anchor_past.reverse_match{
        return f64::MIN;
    }
    if anchor_curr.ref_pos == anchor_past.ref_pos ||
    anchor_curr.query_pos == anchor_past.query_pos{
        return f64::MIN;
    }
    let d_r = (anchor_curr.ref_pos as f64 - anchor_past.ref_pos as f64).abs();
    let d_q;
    if anchor_curr.reverse_match {
        d_q = anchor_past.query_pos as f64 - anchor_curr.query_pos as f64;
    } else {
        d_q = anchor_curr.query_pos as f64 - anchor_past.query_pos as f64;
    }

    if d_q <= 0.{
        return f64::MIN;
    }

    let gap = (d_r - d_q).abs();
    if gap > 3000.{
        return f64::MIN;
    }
    return 100. - gap;
}

pub fn chain_seeds(ref_sketch: &Sketch, query_sketch: &Sketch, k: usize, c_adj: usize) {
    let kmer_seeds_ref = &ref_sketch.kmer_seeds;
    let kmer_seeds_query = &query_sketch.kmer_seeds;
    let kmer_seeds_small;
    let kmer_seeds_large;
    if kmer_seeds_ref.len() > kmer_seeds_query.len() {
        kmer_seeds_small = kmer_seeds_query;
        kmer_seeds_large = kmer_seeds_ref;
    } else {
        kmer_seeds_large = kmer_seeds_query;
        kmer_seeds_small = kmer_seeds_ref;
    }

    let mut anchors = vec![];
    for (canon_kmer, small_pos) in kmer_seeds_small.iter() {
        if kmer_seeds_large.contains_key(canon_kmer) {
            let large_pos = &kmer_seeds_large[canon_kmer];
            for pos1 in small_pos {
                for pos2 in large_pos {
                    anchors.push(Anchor::new(
                        &(pos1.0, pos1.2),
                        &(pos2.0, pos1.2),
                        pos1.1 != pos2.1,
                    ));
                }
            }
        }
    }
    anchors.sort();
    dbg!(kmer_seeds_ref.len(), kmer_seeds_query.len(), anchors.len());

    let mut pointer_vec = vec![0; anchors.len()];
    let mut score_vec = vec![0.; anchors.len()];
    let mut chain_part = partition_vec![];
    for i in 0..anchors.len() {
        chain_part.push(i);
        let anchor_curr = &anchors[i];
        let mut best_score = 0.;
        let mut best_prev_index = i;
        for j in (0..i).rev() {
            let anchor_past = &anchors[j];
            if anchor_curr.ref_pos - anchor_past.ref_pos > 50000 {
                break;
            }
            let new_score = score_anchors(anchor_curr, anchor_past) + score_vec[j];
            if new_score == f64::MIN {
                continue;
            }
            if new_score > best_score {
                best_score = new_score;
                best_prev_index = j;
            }
        }
        score_vec[i] = best_score;
        pointer_vec[i] = best_prev_index;
        if best_prev_index != i {
            chain_part.union(i, best_prev_index);
        }
    }

    let mut intervals = vec![];
    for set in chain_part.all_sets() {
        let mut small_chain = false;
        let mut first_iter = true;
        let mut smallest_id = usize::MAX;
        let mut largest_id = usize::MIN;
        let mut num_anchors = 0;
        for (index, _value) in set {
            if first_iter {
                num_anchors = chain_part.len_of_set(index);
                if num_anchors < 3{
                    small_chain = true;
                    break;
                }
                first_iter = false;
            }
            if index < smallest_id {
                smallest_id = index;
            }
            if index > largest_id {
                largest_id = index;
            }
        }
        let score = score_vec[largest_id];
        if small_chain {
            continue;
        }
        let interval_on_ref = (anchors[smallest_id].ref_pos, anchors[largest_id].ref_pos);
        let endpoint1 = anchors[smallest_id].query_pos;
        let endpoint2 = anchors[largest_id].query_pos;
        let interval_on_query = (
            GnPosition::min(endpoint1, endpoint2),
            GnPosition::max(endpoint1, endpoint2),
        );
        intervals.push((score, interval_on_ref, interval_on_query, num_anchors));
    }
    intervals.sort_by(|x, y| y.partial_cmp(&x).unwrap());
    let mut interval_tree = IntervalTree::new();
    let mut good_intervals = vec![];
    for (i,(_score, int_r, int_q, num_anchors)) in intervals.iter().enumerate() {
        let q_interval = (int_q.0)..(int_q.1);
        if interval_tree.find(&q_interval).count() == 0 {
            interval_tree.insert(q_interval, i);
            good_intervals.push((int_q, int_r, num_anchors));
        }
        else{
            let overlaps = interval_tree.find(&q_interval);
            for ol in overlaps{
                let ol_interval = &intervals[*ol.data()].2;
                if ((ol_interval.1 - int_q.0) as f64) < (ol_interval.1 - ol_interval.0) as f64 * 0.02 ||
               ((int_q.1 - ol_interval.0) as f64) < (ol_interval.1 - ol_interval.0) as f64 * 0.02  {
                    good_intervals.push((int_q, int_r, num_anchors));
                }
            }
        }
    }

    let mut query_cov_length = 0;
    let mut effective_length = 0;
    let mut ref_cov_length = 0;
    let mut num_anchors = 0;
    good_intervals.sort_by(|x,y| x.2.cmp(&y.2));
    for interval in good_intervals{
        let add_query_length = interval.0.1 - interval.0.0 + c_adj as u32;
        let add_ref_length = interval.1.1 - interval.1.0 + c_adj as u32;
        query_cov_length += add_query_length;
        ref_cov_length += add_ref_length;
        effective_length += u32::max(add_query_length, add_ref_length);
        num_anchors += interval.2;
        let ani_est = f64::powf(*interval.2 as f64 * c_adj as f64 / u32::max(add_query_length,add_ref_length) as f64, 1. / k as f64);
//        dbg!(ani_est, interval.2, add_query_length, add_ref_length, interval);

    }

    dbg!(query_cov_length as f64 / query_sketch.total_sequence_length as f64);
    dbg!(ref_cov_length as f64 / ref_sketch.total_sequence_length as f64);

    let ani_est = f64::powf(num_anchors as f64 / effective_length as f64 * c_adj as f64, 1. / k as f64);
    dbg!(ani_est);
}
