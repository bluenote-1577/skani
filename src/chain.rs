use crate::params::*;
use crate::types::*;
use bio::data_structures::interval_tree::IntervalTree;
use fxhash::FxHashMap;
use partitions::*;
use std::mem;

fn wilson_interval(n: f64, n_s: f64, k: usize) -> (f64, f64) {
    let z = 2.0;
    let p_est = n_s / n;
    let shift_mean = 1. / (1. + f64::powi(z, 2) / n) * (p_est + f64::powi(z, 2) / (2. * n));
    let dev_term = z / 1. / (1. + f64::powi(z, 2) / n);
    let dev_term =
        dev_term * f64::sqrt(p_est * (1. - p_est) / n + f64::powi(z, 2) / (4. * f64::powi(n, 2)));
    let upper = f64::powf(shift_mean + dev_term, 1. / k as f64);
    let lower = f64::powf(shift_mean - dev_term, 1. / k as f64);
    let mid = f64::powf(shift_mean, 1. / k as f64);
    return (upper - mid, mid - lower);
}

pub fn map_params_from_sketch(ref_sketch: &Sketch, mode: &str) -> MapParams{
    let fragment_length = fragment_length_formula(ref_sketch.total_sequence_length);
    let max_gap_length = 50.;
    let anchor_score = 50.;
    let min_anchors = 5;
    let length_cutoff = fragment_length;
    let frac_cover_cutoff = 0.001;
    let mode = mode.to_string();
    let chain_band = 200;
    let valid_ks:Vec<(_,_)> = ref_sketch.kmer_seeds_k.iter().enumerate().filter(|x| x.1.len() > 0).collect();
    let k = valid_ks.iter().max_by(|x,y| x.0.cmp(&y.0)).unwrap().0;
    return MapParams{fragment_length, max_gap_length, anchor_score, min_anchors, length_cutoff, mode, frac_cover_cutoff, chain_band, k};

}

pub fn chain_seeds(ref_sketch: &Sketch, query_sketch: &Sketch, map_params: MapParams) {
    let anchor_chunks = get_anchors(ref_sketch, query_sketch, &map_params);
    let chain_results = chain_anchors_ani(&anchor_chunks, &map_params);
    let mut good_intervals_chunks = vec![];
    for i in 0..anchor_chunks.chunks.len() {
        let chain_result = &chain_results[i];
        let anchors = &anchor_chunks.chunks[i];
        let intervals = get_chain_intervals(chain_result, anchors, &map_params);
        let good_intervals = get_nonoverlapping_chains(&intervals, query_sketch.contigs.len());
        good_intervals_chunks.push(good_intervals);
    }
    calculate_ani(
        &good_intervals_chunks,
        ref_sketch,
        query_sketch,
        &anchor_chunks,
        &map_params
    );
}

fn calculate_ani(
    int_chunks: &Vec<Vec<ChainInterval>>,
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
    anchor_chunks: &AnchorChunks,
    map_params: &MapParams
) {
    let k = map_params.k;
    let mut ani_ests = vec![];
    let mut num_good_chunks = 0;
    let mut all_anchors_total = 0;
    let mut total_query_range = 0;
    let mut num_seeds_total = 0;
    for (i, intervals) in int_chunks.iter().enumerate() {
        let mut max_anchors = 0;
        let mut total_anchors = 0;
        let mut total_bases_contained = 0;
        let mut total_range_query = (GnPosition::MAX, GnPosition::MIN);
        let mut total_range_ref = (GnPosition::MAX, GnPosition::MIN);
        for int in intervals {
//                                    dbg!(&int);
            //            if int.num_anchors < max_anchors / 10 {
            //                continue;
            //            } else if int.num_anchors > max_anchors {
            //                max_anchors = int.num_anchors;
            //            }
            total_anchors += int.num_anchors;
            if int.interval_on_query.0 < total_range_query.0 {
                total_range_query.0 = int.interval_on_query.0;
            }
            if int.interval_on_query.1 > total_range_query.1 {
                total_range_query.1 = int.interval_on_query.1;
            }
            if int.interval_on_ref.0 < total_range_ref.0 {
                total_range_ref.0 = int.interval_on_ref.0;
            }
            if int.interval_on_ref.1 > total_range_ref.1 {
                total_range_ref.1 = int.interval_on_ref.1;
            }
            total_bases_contained += int.interval_on_query.1 - int.interval_on_query.0;
        }
//                                dbg!(anchor_chunks.seeds_in_chunk[i]);
        if total_anchors == 0 {
            continue;
        }
        let ml_hits = total_anchors as f64 / anchor_chunks.seeds_in_chunk[i] as f64;
        let ani_est = f64::powf(ml_hits, 1. / k as f64);

        ani_ests.push(ani_est);
//        dbg!(
//            ani_est,
//            total_range_query,
//            total_anchors,
//            anchor_chunks.seeds_in_chunk[i],
//            &intervals 
//        );
        all_anchors_total += total_anchors;
        num_seeds_total += anchor_chunks.seeds_in_chunk[i];
        num_good_chunks += 1;
//        total_query_range += total_range_query.1 - total_range_query.0;
        total_query_range += total_bases_contained;
    }
    let final_ani = ani_ests.iter().sum::<f64>() / ani_ests.len() as f64;

    let (upper, lower) = wilson_interval(num_seeds_total as f64, all_anchors_total as f64, k);
    let covered = total_query_range as f64 / query_sketch.total_sequence_length as f64;
    let q_string;
    if map_params.mode == "classify" {
        q_string = &query_sketch.contigs[0];
    } else {
        q_string = &query_sketch.file_name;
    }
    println!(
        "Query {} Ref {} - ANI {}, +/- = {}/{}. Covered {}",
        q_string,
        ref_sketch.file_name,
        final_ani,
        final_ani + upper,
        final_ani - lower,
        covered,
    );

    //    let debias_2 =  eff_len * f64::powf(final_ani, k as f64) * term2;
    //    let debias_2_final = (final_ani - f64::powf(1./eff_len, 1./k as f64)*debias_2) * (debias - debias_2);
}

pub fn score_anchors(anchor_curr: &Anchor, anchor_past: &Anchor, map_params: &MapParams) -> f64 {
    if anchor_curr.ref_contig != anchor_past.ref_contig {
        return f64::MIN;
    }
    if anchor_curr.reverse_match != anchor_past.reverse_match {
        return f64::MIN;
    }

    if anchor_curr.ref_pos == anchor_past.ref_pos || anchor_curr.query_pos == anchor_past.query_pos
    {
        return f64::MIN;
    }
    let d_q = (anchor_curr.query_pos as f64 - anchor_past.query_pos as f64).abs();
    let d_r;
    if anchor_curr.reverse_match {
        d_r = anchor_past.ref_pos as f64 - anchor_curr.ref_pos as f64;
    } else {
        d_r = anchor_curr.ref_pos as f64 - anchor_past.ref_pos as f64;
    }

    if d_r <= 0. {
        return f64::MIN;
    }

    let gap = (d_r - d_q).abs();
    if gap > map_params.max_gap_length{
        return f64::MIN;
    }
    return map_params.anchor_score - gap;
}

fn get_anchors(ref_sketch: &Sketch, query_sketch: &Sketch, map_params: &MapParams) -> AnchorChunks {
    let k = map_params.k;
    let kmer_seeds_ref = &ref_sketch.kmer_seeds_k[k];
    let kmer_seeds_query = &query_sketch.kmer_seeds_k[k];
    let mut anchors = vec![];
    let mut count_vec = vec![];
    for ref_pos in kmer_seeds_ref.values(){
        count_vec.push(ref_pos.len());
    }
    count_vec.sort();
    let mut max_repet_cutoff = count_vec[count_vec.len() - count_vec.len() / 10000 - 1];
    if max_repet_cutoff< 30{
        max_repet_cutoff = usize::MAX;
    }
    dbg!(max_repet_cutoff);
    let mut query_kmers_with_hits = 0;
    let mut query_positions_all = vec![vec![]; query_sketch.contigs.len()];
    for (canon_kmer, query_pos) in kmer_seeds_query.iter() {
//        dbg!(KmerEnc::print_string(canon_kmer.kmer,k));
        if query_pos.len() > max_repet_cutoff {
            continue;
        }

        for qpos in query_pos.iter() {
            query_positions_all[qpos.2 as usize].push(qpos.0);
        }
        if kmer_seeds_ref.contains_key(canon_kmer) {
            let ref_pos = &kmer_seeds_ref[canon_kmer];

            if ref_pos.len() > max_repet_cutoff{
                continue;
            }
            query_kmers_with_hits += 1;
            for qpos in query_pos {
                for rpos in ref_pos {
                    anchors.push(Anchor::new(
                        &(rpos.0, rpos.2),
                        &(qpos.0, qpos.2),
                        rpos.1 != qpos.1,
                    ));
                }
            }
        }
        else{
//            KmerEnc::print_string(canon_kmer.kmer, k);
//            dbg!(query_pos);
        }
    }
    if anchors.is_empty() {
        return AnchorChunks::default();
    }
    anchors.sort();
    for query_position_vec in query_positions_all.iter_mut() {
        query_position_vec.sort();
    }
    dbg!(
        kmer_seeds_ref.len(),
        kmer_seeds_query.len(),
        anchors.len(),
        query_kmers_with_hits,
        f64::powf((query_kmers_with_hits as f64)/ (kmer_seeds_query.len() as f64), 1./(k as f64)),
    );
    //        dbg!(&anchors);
    let mut lengths = vec![];
    let mut chunks = vec![];
    let mut curr_anchor_chunk = vec![];
    let mut block_seeds = vec![];
    let smallest_anchor_query_pos = anchors[0].query_pos;
    let mut last_query_contig = anchors[0].query_contig;
    let mut curr_end_point = smallest_anchor_query_pos + map_params.fragment_length as u32;
    let mut running_counter = 0;
    for anchor in anchors {
        if last_query_contig != anchor.query_contig
            || anchor.query_pos > curr_end_point
        {
            if query_positions_all[last_query_contig as usize].is_empty() {
                dbg!(&query_sketch.contigs[last_query_contig as usize]);
                continue;
            }
            let mut num_seeds_in_block = 1;
            loop {
                if running_counter >= query_positions_all[last_query_contig as usize].len() {
                    break;
                }
                if query_positions_all[last_query_contig as usize][running_counter]
                    < curr_end_point - map_params.fragment_length as u32
                {
                    running_counter += 1;
                } else if (query_positions_all[last_query_contig as usize][running_counter]
                    <= curr_end_point)
                    || (last_query_contig != anchor.query_contig)
                {
                    running_counter += 1;
                    num_seeds_in_block += 1;
                } else {
                    break;
                }
            }
            block_seeds.push(num_seeds_in_block);
            curr_end_point += map_params.fragment_length as u32;
            chunks.push(mem::take(&mut curr_anchor_chunk));
            curr_anchor_chunk = vec![];
            lengths.push(map_params.fragment_length as u32);
            if last_query_contig != anchor.query_contig {
                curr_end_point = anchor.query_pos + map_params.fragment_length as u32;
                running_counter = 0;
            }
        }
        last_query_contig = anchor.query_contig;
        curr_anchor_chunk.push(anchor);
    }
    if curr_anchor_chunk.len() > 0 {
        let mut num_seeds_in_block = 1;
        loop {
            if query_positions_all[last_query_contig as usize].is_empty() {
                dbg!(&query_sketch.contigs[last_query_contig as usize]);
                continue;
            }
            if running_counter >= query_positions_all[last_query_contig as usize].len() {
                break;
            }
            if query_positions_all[last_query_contig as usize][running_counter]
                < curr_end_point - map_params.fragment_length as u32
            {
                running_counter += 1;
            } else if (query_positions_all[last_query_contig as usize][running_counter]
                <= curr_anchor_chunk[curr_anchor_chunk.len() - 1].query_pos)
                || (last_query_contig
                    != curr_anchor_chunk[curr_anchor_chunk.len() - 1].query_contig)
            {
                running_counter += 1;
                num_seeds_in_block += 1;
            } else {
                break;
            }
        }
        lengths.push(
            curr_anchor_chunk[curr_anchor_chunk.len() - 1].query_pos
                - curr_anchor_chunk[0].query_pos,
        );
        chunks.push(mem::take(&mut curr_anchor_chunk));
        block_seeds.push(num_seeds_in_block);
    }

    return AnchorChunks {
        chunks,
        lengths,
        seeds_in_chunk: block_seeds,
    };
}

fn chain_anchors_ani(anchor_chunks: &AnchorChunks, map_params: &MapParams) -> Vec<ChainingResult> {
    let mut chaining_results = vec![];

    let num_chunks = anchor_chunks.chunks.len();
    let past_chain_length = usize::min(map_params.fragment_length / 2, 5000);
    for anchor_chunk in anchor_chunks.chunks.iter() {
        let mut pointer_vec = vec![0; anchor_chunk.len()];
        let mut score_vec = vec![0.; anchor_chunk.len()];
        let mut chain_part = partition_vec![];

        for i in 0..anchor_chunk.len() {
            chain_part.push(i);
            let anchor_curr = &anchor_chunk[i];
            let mut best_score = 0.;
            let mut best_prev_index = i;
            for j in (0..i).rev() {
                let anchor_past = &anchor_chunk[j];
                if anchor_curr.query_contig != anchor_past.query_contig {
                    continue;
                }
                if anchor_curr.ref_contig != anchor_past.ref_contig {
                    continue;
                }
                if anchor_curr.query_pos - anchor_past.query_pos > past_chain_length as u32 ||
                i - j > map_params.chain_band{
                    break;
                }
                let anchor_score = score_anchors(anchor_curr, anchor_past, map_params);
                if anchor_score == f64::MIN {
                    continue;
                }
                let new_score = anchor_score + score_vec[j];
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

        chaining_results.push(ChainingResult {
            pointer_vec,
            chain_part,
            score_vec,
            num_chunks,
        });
    }
    return chaining_results;
}

//fn chain_anchors_global(anchors: &Vec<Anchor>) -> ChainingResult {
//    let mut pointer_vec = vec![0; anchors.len()];
//    let mut score_vec = vec![0.; anchors.len()];
//    let mut chain_part = partition_vec![];
//    for i in 0..anchors.len() {
//        chain_part.push(i);
//        let anchor_curr = &anchors[i];
//        let mut best_score = 0.;
//        let mut best_prev_index = i;
//        for j in (0..i).rev() {
//            let anchor_past = &anchors[j];
//            if anchor_curr.query_contig != anchor_past.query_contig {
//                break;
//            }
//            if anchor_curr.query_pos - anchor_past.query_pos > u32::min(FRAGMENT_LENGTH as u32, 2000) {
//                break;
//            }
//            let new_score = score_anchors(anchor_curr, anchor_past) + score_vec[j];
//            if new_score == f64::MIN {
//                continue;
//            }
//            if new_score > best_score {
//                best_score = new_score;
//                best_prev_index = j;
//            }
//        }
//        score_vec[i] = best_score;
//        pointer_vec[i] = best_prev_index;
//        if best_prev_index != i {
//            chain_part.union(i, best_prev_index);
//        }
//    }
//
//    return ChainingResult {
//        pointer_vec,
//        chain_part,
//        score_vec,
//        num_chunks: usize::MAX,
//    };
//}

fn get_chain_intervals(cr: &ChainingResult, anchors: &Vec<Anchor>, map_params: &MapParams) -> Vec<ChainInterval> {
    let mut intervals = vec![];
    for set in cr.chain_part.all_sets() {
        let mut small_chain = false;
        let mut first_iter = true;
        let mut smallest_id = usize::MAX;
        let mut largest_id = usize::MIN;
        let mut num_anchors = 0;
        for (index, _value) in set {
            if first_iter {
                num_anchors = cr.chain_part.len_of_set(index);
                if num_anchors < map_params.min_anchors {
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
            assert!((anchors[smallest_id].ref_contig == anchors[index].ref_contig));
            assert!((anchors[smallest_id].query_contig == anchors[index].query_contig));
        }
        let score = cr.score_vec[largest_id];
        if small_chain {
            continue;
        }
        let interval_on_query = (
            anchors[smallest_id].query_pos,
            anchors[largest_id].query_pos,
        );
        let endpoint1 = anchors[smallest_id].ref_pos;
        let endpoint2 = anchors[largest_id].ref_pos;
        let interval_on_ref = (
            GnPosition::min(endpoint1, endpoint2),
            GnPosition::max(endpoint1, endpoint2),
        );
        let ref_contig = anchors[smallest_id].ref_contig as usize;
        let query_contig = anchors[smallest_id].query_contig as usize;
        let chain_interval = ChainInterval {
            interval_on_query,
            interval_on_ref,
            ref_contig,
            query_contig,
            score,
            num_anchors,
        };
        intervals.push(chain_interval);
    }
    intervals.sort_by(|y, x| x.partial_cmp(&y).unwrap());
    intervals
}
fn get_nonoverlapping_chains(
    intervals: &Vec<ChainInterval>,
    num_contigs: usize,
) -> Vec<ChainInterval> {
    let mut interval_trees = FxHashMap::default();
    let mut good_intervals = vec![];
    for (i, int) in intervals.iter().enumerate() {
        let q_interval = (int.interval_on_query.0)..(int.interval_on_query.1);
        let contig_id = int.query_contig;
        let interval_tree = interval_trees
            .entry(contig_id)
            .or_insert(IntervalTree::new());
        if interval_tree.find(&q_interval).count() == 0 {
            interval_tree.insert(q_interval, i);
            good_intervals.push(int.clone());
        } else {
            let overlaps = interval_tree.find(&q_interval);
            let mut small_ol = true;
            for ol in overlaps {
                let ol_interval = &intervals[*ol.data()];
                let interval_length = ol_interval.query_range_len();
                let overlap = GnPosition::min(
                    int.interval_on_query.1 - ol_interval.interval_on_query.0,
                    ol_interval.interval_on_query.1 - int.interval_on_query.0,
                );
                if (overlap as f64) > -1. {
                    //interval_length as f64 * 0.00 {
                    small_ol = false;
                }
            }
            if small_ol {
                good_intervals.push(int.clone());
            }
        }
    }
    good_intervals.sort_by(|x, y| y.num_anchors.cmp(&x.num_anchors));
    return good_intervals;
}

fn calculate_ani_chain(
    good_intervals: &Vec<&ChainInterval>,
    c_adj: f64,
    k: usize,
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
) {
    let mut query_cov_length = 0;
    let mut effective_length = 0;
    let mut ref_cov_length = 0;
    let mut num_anchors = 0;
    //    good_intervals.sort_by(|x, y| x.3.cmp(&y.3));
    for interval in good_intervals {
        let add_query_length = interval.query_range_len() + 2 * c_adj as u32;
        let add_ref_length = interval.ref_range_len() + 2 * c_adj as u32;
        query_cov_length += add_query_length;
        ref_cov_length += add_ref_length;
        effective_length += u32::max(add_query_length, add_ref_length);
        num_anchors += interval.num_anchors;
        let ani_est = f64::powf(
            interval.num_anchors as f64 * c_adj as f64
                / u32::max(add_query_length, add_ref_length) as f64,
            1. / k as f64,
        );
        //        dbg!(
        //            ani_est,
        //            interval.num_anchors,
        //            add_query_length,
        //            add_ref_length,
        //            interval
        //        );
    }

    dbg!(query_cov_length as f64 / query_sketch.total_sequence_length as f64);
    dbg!(ref_cov_length as f64 / ref_sketch.total_sequence_length as f64);

    let ani_est = f64::powf(
        num_anchors as f64 / effective_length as f64 * c_adj as f64,
        1. / k as f64,
    );
    dbg!(ani_est);
}
