use crate::params::*;
use crate::types::*;
use bio::data_structures::interval_tree::IntervalTree;
use fxhash::FxHashMap;
use log::{debug, info, trace, warn};
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

pub fn map_params_from_sketch(
    ref_sketch: &Sketch,
    mode: &str,
    amino_acid: bool,
    euk: bool,
    screen: bool,
) -> MapParams {
    let fragment_length =
        fragment_length_formula(ref_sketch.total_sequence_length, amino_acid, euk);
    let max_gap_length = D_MAX_GAP_LENGTH;
    let anchor_score = 50.;
    let min_anchors = 4;
    let length_cutoff = fragment_length;
    trace!("Fragment length is {}.", fragment_length);
    let frac_cover_cutoff;
    if amino_acid {
        frac_cover_cutoff = D_FRAC_COVER_CUTOFF_AA;
    } else {
        frac_cover_cutoff = D_FRAC_COVER_CUTOFF;
    }
    let length_cover_cutoff = 5000000;
    let mode = mode.to_string();
    let chain_band = D_CHAIN_BAND;
    let min_score = min_anchors as f64 * anchor_score;
    let valid_ks: Vec<(_, _)> = ref_sketch
        .kmer_seeds_k
        .iter()
        .enumerate()
        .filter(|x| x.1.len() > 0)
        .collect();
    let k = valid_ks.iter().max_by(|x, y| x.0.cmp(&y.0)).unwrap().0;
    return MapParams {
        fragment_length,
        max_gap_length,
        anchor_score,
        min_anchors,
        length_cutoff,
        mode,
        frac_cover_cutoff,
        length_cover_cutoff,
        chain_band,
        k,
        amino_acid,
        min_score,
        euk,
        screen,
    };
}

pub fn chain_seeds(
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
    map_params: MapParams,
) -> AniEstResult {
    let (anchor_chunks, switched) = get_anchors(ref_sketch, query_sketch, &map_params);
    let chain_results = chain_anchors_ani(&anchor_chunks, &map_params);
    let mut good_intervals = vec![];
    for i in 0..anchor_chunks.chunks.len() {
        let chain_result = &chain_results[i];
        let anchors = &anchor_chunks.chunks[i];
        get_chain_intervals(&mut good_intervals, chain_result, anchors, &map_params, i);
    }
    let good_interval_chunks =
        get_nonoverlapping_chains(&mut good_intervals, anchor_chunks.chunks.len());
    return calculate_ani(
        &good_interval_chunks,
        ref_sketch,
        query_sketch,
        &anchor_chunks,
        &map_params,
        switched,
    );
}

fn calculate_ani(
    int_chunks: &Vec<Vec<ChainInterval>>,
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
    anchor_chunks: &AnchorChunks,
    map_params: &MapParams,
    switched: bool,
) -> AniEstResult {
    let k = map_params.k;
    let mut ani_ests = vec![];
    let mut num_good_chunks = 0;
    let mut all_anchors_total = 0;
    let mut total_query_range = 0;
    let mut total_ref_range = 0;
    let mut num_seeds_total = 0;
    for (i, intervals) in int_chunks.iter().enumerate() {
        let mut max_anchors = 0;
        let mut total_anchors = 0;
        let mut total_bases_contained_query = 0;
        let mut total_bases_contained_ref = 0;
        let mut total_range_query = (GnPosition::MAX, GnPosition::MIN);
        let mut total_range_ref = (GnPosition::MAX, GnPosition::MIN);
        for int in intervals {
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
            if !switched {
                total_bases_contained_query += int.interval_on_query.1 - int.interval_on_query.0
                    + map_params.k as GnPosition
                    + 2 * ref_sketch.c as GnPosition;
                total_bases_contained_ref += int.interval_on_ref.1 - int.interval_on_ref.0
                    + map_params.k as GnPosition
                    + 2 * ref_sketch.c as GnPosition;
            } else {
                total_bases_contained_query += int.interval_on_ref.1 - int.interval_on_ref.0
                    + map_params.k as GnPosition
                    + 2 * ref_sketch.c as GnPosition;
                total_bases_contained_ref += int.interval_on_query.1 - int.interval_on_query.0
                    + map_params.k as GnPosition
                    + 2 * ref_sketch.c as GnPosition;
            }
        }

        if total_anchors == 0 {
            continue;
        }
        let ml_hits = if map_params.amino_acid {
            f64::min(
                1.,
                total_anchors as f64 / anchor_chunks.seeds_in_chunk[i] as f64 * 6.,
            )
        } else {
            f64::min(
                1.,
                total_anchors as f64 / anchor_chunks.seeds_in_chunk[i] as f64,
            )
        };
        let ani_est = if map_params.amino_acid {
            f64::powf(ml_hits, 1. / k as f64)
        } else {
            f64::powf(ml_hits, 1. / k as f64)
        };
        total_query_range += total_bases_contained_query + 2 * ref_sketch.c as GnPosition;
        total_ref_range += total_bases_contained_ref + 2 * ref_sketch.c as GnPosition;

        //        total_bases_contained_query =
        //            total_range_query.1 - total_range_query.0 + map_params.k as GnPosition;
        //        total_bases_contained_ref =
        //            total_range_query.1 - total_range_query.0 + map_params.k as GnPosition;
        //
        //        total_query_range += total_bases_contained_query;
        //        total_ref_range += total_bases_contained_ref;

        //        ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i]));
        //        ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i]));
        ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i]));
        trace!(
            "Ani est fragment {}, total range {:?}, total anchors {}, seeds in fragment {:?}",
            ani_est,
            total_range_query,
            total_anchors,
            anchor_chunks.seeds_in_chunk[i]
        );
        trace!(
            "Intervals {:?}, Num Anchors in Interval {}",
            &intervals,
            intervals[0].num_anchors
        );
        all_anchors_total += total_anchors;
        num_seeds_total += anchor_chunks.seeds_in_chunk[i];
        num_good_chunks += 1;
    }
    ani_ests.sort_by(|x, y| x.partial_cmp(&y).unwrap());

//    for ani in ani_ests.iter(){
//        println!("{},{}", ani.0, ani.1);
//    }

    if ani_ests.len() == 0 {
        let mut ret = AniEstResult::default();
        ret.ani = f64::NAN;
        return ret;
    }
    let lower = 10 * ani_ests.len() / 100;
    let upper = 90 * ani_ests.len() / 100;
    let mut weighted_avg = 0.;
    let mut total_weight_interval = 0;
    for i in lower..upper + 1 {
        weighted_avg += ani_ests[i].0 * (ani_ests[i].1 as f64);
        total_weight_interval += ani_ests[i].1;
    }
    let mut final_ani = weighted_avg / total_weight_interval as f64;
    //    let final_ani = quantile_slice.iter().sum::<f64>() / quantile_slice.len() as f64;

    let (upper, lower) = wilson_interval(num_seeds_total as f64, all_anchors_total as f64, k);
    let covered_query = f64::min(1., total_query_range as f64 / query_sketch.total_sequence_length as f64);
    let covered_ref = f64::min(1.,total_ref_range as f64 / ref_sketch.total_sequence_length as f64);
    let q_string;
    if map_params.mode == "classify" {
        q_string = &query_sketch.contigs[0];
    } else {
        q_string = &query_sketch.file_name;
    }
    let id_string = if map_params.amino_acid { "AAI" } else { "ANI" };
    debug!(
        "Query {} Ref {} - {} {}, +/- = {}/{}. Covered {}",
        q_string,
        ref_sketch.file_name,
        id_string,
        final_ani,
        final_ani + upper,
        final_ani - lower,
        covered_query,
    );

    if covered_query < map_params.frac_cover_cutoff
        && total_query_range < map_params.length_cover_cutoff as GnPosition
    {
        final_ani = -1.;
    }
    return AniEstResult {
        ani: final_ani,
        align_fraction_query: covered_query,
        align_fraction_ref: covered_ref,
        ref_file: ref_sketch.file_name.clone(),
        query_file: query_sketch.file_name.clone(),
        query_contig: query_sketch.contigs[0].clone(),
    };
}

pub fn score_anchors(anchor_curr: &Anchor, anchor_past: &Anchor, map_params: &MapParams) -> f64 {
    if anchor_curr.query_phase != anchor_past.query_phase
        || anchor_curr.ref_phase != anchor_past.ref_phase
    {
        return f64::MIN;
    }
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
    let acqpf64 = anchor_curr.query_pos as f64;
    let apqpf64 = anchor_past.query_pos as f64;
    let acrpf64 = anchor_curr.ref_pos as f64;
    let aprpf64 = anchor_past.ref_pos as f64;

    let d_q = (acqpf64 - apqpf64).abs();
    let d_r;
    if anchor_curr.reverse_match {
        d_r = aprpf64 - acrpf64;
    } else {
        d_r = acrpf64 - aprpf64;
    }

    if d_r <= 0. {
        return f64::MIN;
    }

    let gap = (d_r - d_q).abs();
    if gap > map_params.max_gap_length {
        return f64::MIN;
    }
//    let ol_q = f64::max(0., map_params.k as f64 - d_q);
//    let ol_r = f64::max(0., map_params.k as f64 - d_r);
//    let ol = f64::max(ol_q, ol_r);
//    return map_params.anchor_score * (1. - ol / map_params.k as f64) - gap;
    return  map_params.anchor_score - gap
}

fn get_anchors(
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
    map_params: &MapParams,
) -> (AnchorChunks, bool) {
    let k = map_params.k;
    let kmer_seeds_ref;
    let kmer_seeds_query;
    let mut query_positions_all;
    let switched;
    //    if query_sketch.contigs.len() < ref_sketch.contigs.len(){
    if query_sketch.total_sequence_length > ref_sketch.total_sequence_length {
        switched = true;
        kmer_seeds_ref = &query_sketch.kmer_seeds_k[k];
        kmer_seeds_query = &ref_sketch.kmer_seeds_k[k];
        query_positions_all = vec![vec![]; ref_sketch.contigs.len()];
    } else {
        switched = false;
        kmer_seeds_ref = &ref_sketch.kmer_seeds_k[k];
        kmer_seeds_query = &query_sketch.kmer_seeds_k[k];
        query_positions_all = vec![vec![]; query_sketch.contigs.len()];
    }
    //    let kmer_seeds_ref = &ref_sketch.kmer_seeds_k[k];
    //    let kmer_seeds_query = &query_sketch.kmer_seeds_k[k];
    let mut anchors = vec![];
    let mut query_kmers_with_hits = 0;
    for (canon_kmer, query_pos) in kmer_seeds_query.iter() {
        //        dbg!(KmerEnc::print_string(canon_kmer.kmer,k));
        if query_pos.len() > query_sketch.repetitive_kmers[k] {
            continue;
        }
        let contains = kmer_seeds_ref.contains_key(canon_kmer);

        if !contains{
            for qpos in query_pos.iter() {
                query_positions_all[qpos.contig_index as usize].push(qpos.pos);
            }
        }

        else {
            let ref_pos = &kmer_seeds_ref[canon_kmer];

            if ref_pos.len() > ref_sketch.repetitive_kmers[k] {
                continue;
            }

            for qpos in query_pos.iter() {
                query_positions_all[qpos.contig_index as usize].push(qpos.pos);
            }

            query_kmers_with_hits += 1;
            for qpos in query_pos {
                for rpos in ref_pos {
                    anchors.push(Anchor::new(
                        &(rpos.pos, rpos.contig_index),
                        &(qpos.pos, qpos.contig_index),
                        rpos.phase,
                        qpos.phase,
                        rpos.canonical != qpos.canonical,
                    ));
                }
            }
        } 
    }
    if anchors.is_empty() {
        info!(
            "no anchors found for {}, {}",
            &ref_sketch.file_name, &query_sketch.file_name
        );
        return (AnchorChunks::default(), true);
    }
    anchors.sort();
    for query_position_vec in query_positions_all.iter_mut() {
        query_position_vec.sort();
    }
    trace!(
        "Ref seeds len {}, Query seeds len {}, Anchors {}, Seeds hit query {}, Est {}",
        kmer_seeds_ref.len(),
        kmer_seeds_query.len(),
        anchors.len(),
        query_kmers_with_hits,
        f64::powf(
            (query_kmers_with_hits as f64) / (kmer_seeds_query.len() as f64),
            1. / (k as f64)
        ),
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
        if last_query_contig != anchor.query_contig || anchor.query_pos > curr_end_point {
            if query_positions_all[last_query_contig as usize].is_empty() {
                warn!("{}", &query_sketch.contigs[last_query_contig as usize]);
                continue;
            }
            let mut num_seeds_in_block = 0;
            let mut first_iter = true;
            loop {
                if running_counter >= query_positions_all[last_query_contig as usize].len() {
                    break;
                }
                if query_positions_all[last_query_contig as usize][running_counter]
                    <= curr_end_point
                {
                    if first_iter {
                        first_iter = false;
                        trace!("start {}", query_positions_all[last_query_contig as usize][running_counter]);
                    }
                    running_counter += 1;
                    num_seeds_in_block += 1
                } else {
                    trace!("end {}", query_positions_all[last_query_contig as usize][running_counter-1]);
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
        let mut num_seeds_in_block = 0;
        loop {
            if query_positions_all[last_query_contig as usize].is_empty() {
                warn!("{}", &query_sketch.contigs[last_query_contig as usize]);
                continue;
            }
            if running_counter >= query_positions_all[last_query_contig as usize].len() {
                break;
            }
            if (query_positions_all[last_query_contig as usize][running_counter]
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

    assert!(block_seeds.len() == chunks.len());

    return (
        AnchorChunks {
            chunks,
            lengths,
            seeds_in_chunk: block_seeds,
        },
        switched,
    );
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
                if anchor_curr.query_pos - anchor_past.query_pos > past_chain_length as u32
                    || i - j > map_params.chain_band
                {
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

fn get_chain_intervals(
    good_intervals: &mut Vec<ChainInterval>,
    cr: &ChainingResult,
    anchors: &Vec<Anchor>,
    map_params: &MapParams,
    chunk_id: usize,
) {
    for set in cr.chain_part.all_sets() {
        let mut small_chain = false;
        let mut first_iter = true;
        let mut max_score = f64::MIN;
        let mut best_index = usize::MAX;
        let mut num_anchors = 1;
        for (index, _value) in set {
            if first_iter {
                if cr.chain_part.len_of_set(index) < map_params.min_anchors {
                    small_chain = true;
                    break;
                }
                first_iter = false;
            }
            if cr.score_vec[index] > max_score {
                max_score = cr.score_vec[index];
                best_index = index;
            }
        }
        if small_chain {
            continue;
        }

        let mut index = best_index;
        while cr.pointer_vec[index] != index {
            index = cr.pointer_vec[index];
            num_anchors += 1;
        }
        small_chain = num_anchors < map_params.min_anchors;
        if small_chain || max_score < map_params.min_score {
            continue;
        }
        let smallest_id = index;
        let largest_id = best_index;
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
            score: max_score,
            num_anchors,
            chunk_id,
        };
        good_intervals.push(chain_interval);
    }
}
fn get_nonoverlapping_chains(
    intervals: &mut Vec<ChainInterval>,
    num_chunks: usize,
) -> Vec<Vec<ChainInterval>> {
    intervals.sort_by(|x, y| y.partial_cmp(&x).unwrap());
    let mut interval_trees = FxHashMap::default();
    let mut interval_trees_ref = FxHashMap::default();
    let mut good_non_overlap_intervals = vec![vec![]; num_chunks];
    for (i, int) in intervals.iter().enumerate() {
        let q_interval = (int.interval_on_query.0)..(int.interval_on_query.1);
        let r_interval = (int.interval_on_ref.0)..(int.interval_on_ref.1);
        let interval_tree_ref = interval_trees_ref
            .entry(int.ref_contig)
            .or_insert(IntervalTree::new());
        let interval_tree_q = interval_trees
            .entry(int.query_contig)
            .or_insert(IntervalTree::new());

        let no_overlap_ref;
        if interval_tree_ref.find(&r_interval).count() == 0 {
            no_overlap_ref = true
        } else {
            no_overlap_ref = false;
        }

        let no_overlap_query;
        if interval_tree_q.find(&q_interval).count() == 0 {
            no_overlap_query = true;
        } else {
            let overlaps = interval_tree_q.find(&q_interval);
            let mut small_ol = false;
            //            for ol in overlaps {
            //                let ol_interval = &intervals[*ol.data()];
            //                let interval_length = ol_interval.query_range_len();
            //                let overlap = GnPosition::min(
            //                    int.interval_on_query.1 - ol_interval.interval_on_query.0,
            //                    ol_interval.interval_on_query.1 - int.interval_on_query.0,
            //                );
            //                if (overlap as f64) < -1. {
            //                    //interval_length as f64 * 0.00 {
            //                    small_ol = true;
            //                }
            //            }
            if small_ol {
                no_overlap_query = true;
            } else {
                no_overlap_query = false;
            }
        }
        if no_overlap_ref && no_overlap_query {
            //                    if  no_overlap_query{
            //

            interval_tree_q.insert(q_interval, i);
            interval_tree_ref.insert(r_interval, i);
            good_non_overlap_intervals[int.chunk_id].push(int.clone());
        }
    }
    //good_non_overlap_intervals.sort_by(|x, y| y.partial_cmp(&x).unwrap());
    return good_non_overlap_intervals;
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
