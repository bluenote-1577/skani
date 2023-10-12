use crate::params::*;
use gbdt::gradient_boost::GBDT;
use crate::types::*;
use bio::data_structures::interval_tree::IntervalTree;
use crate::regression;

use fxhash::FxHashMap;
use log::*;
use partitions::*;
use std::mem;
extern crate interval;
use gcollections::ops::set::*;
use interval::interval_set::*;

fn switch_qr(med_ctg_len_r: f64, med_ctg_len_q: f64, q_sk_len: f64,r_sk_len: f64, query_file_name: &str, ref_file_name: &str)-> bool{
    let score_query = q_sk_len
        * (f64::min(med_ctg_len_q, 300000.));
    let score_ref = r_sk_len
        * (f64::min(med_ctg_len_r, 300000.));
    if score_query == score_ref{
        query_file_name > ref_file_name
    }
    else{
        score_query > score_ref
    }
}

fn mean(data: &[f64]) -> Option<f64> {
    let sum = data.iter().sum::<f64>() as f64;
    let count = data.len() as f64;
    if data.len() > 0{
        return Some(sum/count);
    }
    else{
        return None;
    }
}

fn std_deviation(data: &[f64]) -> f64 {
    let count = data.len();
    let data_mean = mean(&data);
    if data_mean.is_none(){
        return 0.
    }
    else{
        let data_mean = data_mean.unwrap();
        let variance = data.iter().map(|value| {
            let diff = data_mean - (*value);

            diff * diff
        }).sum::<f64>() / count as f64;

        variance.sqrt()
    }
}

fn bootstrap_interval(ani_ests: &Vec<(f64,usize)>) -> (f64,f64,f64){
    let ani_est_no_mult = ani_ests.iter().map(|x| x.0).collect::<Vec<f64>>();
    let std = std_deviation(&ani_est_no_mult);
    let mut res = vec![];
    let mut mult_ani_ests = vec![];
    fastrand::seed(7);
    let num_samp = ani_ests.len();
    //Return no confidence interval if number of samples is too small. 
    if num_samp < 10 {
        return (0.,1., std);
    }
    for (ani,mult) in ani_ests.iter(){
        for _ in 0..*mult{
            mult_ani_ests.push(ani);
        }
    }
    let iters = 100;
    for _ in 0..iters{
        let mut rand_vec = vec![];
        rand_vec.reserve(num_samp);
        for _ in 0..num_samp{
            rand_vec.push(fastrand::usize(..mult_ani_ests.len()));
        }
        let sum = rand_vec.into_iter().map(|x| mult_ani_ests[x]).sum::<f64>();
        res.push(sum/(num_samp as f64));
    }
    res.sort_by(|x,y| x.partial_cmp(y).unwrap());
    (res[iters * 5 / 100 - 1],res[iters * 95 / 100 - 1], std)

}

pub fn map_params_from_sketch <'a>(
    ref_sketch: &Sketch,
    amino_acid: bool,
    command_params: &CommandParams,
    model_opt: &'a Option<GBDT>
) -> MapParams<'a> {
    let max_gap_length = if amino_acid{D_MAX_GAP_LENGTH_AAI} else {D_MAX_GAP_LENGTH};
    let anchor_score = if amino_acid{D_ANCHOR_SCORE_AAI} else {D_ANCHOR_SCORE_ANI};
    let min_anchors = if amino_acid{D_MIN_ANCHORS_AAI} else {D_MIN_ANCHORS_ANI};
    let min_length_cover = if amino_acid{MIN_LENGTH_COVER_AAI} else {MIN_LENGTH_COVER};
    let fragment_length = fragment_length_formula(ref_sketch.total_sequence_length, amino_acid);
    let length_cutoff = fragment_length;
    let mut frac_cover_cutoff = command_params.min_aligned_frac;
    if frac_cover_cutoff < 0.{
        if amino_acid {
            frac_cover_cutoff = D_FRAC_COVER_CUTOFF_AA.parse::<f64>().unwrap()/100.;
        } else {
            frac_cover_cutoff = D_FRAC_COVER_CUTOFF.parse::<f64>().unwrap()/100.;
        }
    }
    let length_cover_cutoff = 5000000;
    let bp_chain_band = if amino_acid {BP_CHAIN_BAND_AAI} else {BP_CHAIN_BAND};
    let index_chain_band = bp_chain_band/ref_sketch.c;
    let min_score = min_anchors as f64 * anchor_score * 0.75;
//    let min_score = 0.;
    let k = ref_sketch.k;
    let model;
    if let Some(m) = model_opt{
        model = Some(m);
    }
    else{
        model = None
    }
    MapParams {
        fragment_length,
        max_gap_length,
        anchor_score,
        min_anchors,
        length_cutoff,
        frac_cover_cutoff,
        length_cover_cutoff,
        index_chain_band,
        k,
        amino_acid,
        min_score,
        robust: command_params.robust,
        median: command_params.median,
        bp_chain_band,
        min_length_cover,
        model
    }
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
    let mut ani = calculate_ani(
        &good_interval_chunks,
        ref_sketch,
        query_sketch,
        &anchor_chunks,
        &map_params,
        switched,
    );
    if let Some(model) = map_params.model{
        regression::predict_from_ani_res(&mut ani, model);
    }
    ani
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
    let c = ref_sketch.c as GnPosition;
    let sensitive_af;
    if c < 200{
        sensitive_af = true;
    }
    else{
        sensitive_af = false;
    }
    let mut _num_good_chunks = 0;
    let mut _all_anchors_total = 0;
    let mut total_query_bases = 0;
    let mut total_ref_range = 0;
    let mut leftmost_interval = &ChainInterval::default();
    let mut rightmost_interval = &ChainInterval::default();
    let mut avg_chain_int_len = 0;
    let mut num_chains = 0;
    for (i, intervals) in int_chunks.iter().enumerate() {
        let mut all_intervals = vec![].to_interval_set();
        let mut total_anchors = 0;
        let mut total_bases_contained_query = 0;
        let mut _total_bases_contained_ref = 0;
        let mut total_range_query = (GnPosition::MAX, GnPosition::MIN);
        let mut total_range_ref = (GnPosition::MAX, GnPosition::MIN);
        for int in intervals {
            total_anchors += int.num_anchors;

            if int.interval_on_query.0 < total_range_query.0 {
                total_range_query.0 = int.interval_on_query.0;
                leftmost_interval = int;
            }
            if int.interval_on_query.1 > total_range_query.1 {
                total_range_query.1 = int.interval_on_query.1;
                rightmost_interval = int;
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
                    + 2 * c;
                _total_bases_contained_ref += int.interval_on_ref.1 - int.interval_on_ref.0
                    + map_params.k as GnPosition
                    + 2 * c;
            } else {
                total_bases_contained_query += int.interval_on_ref.1 - int.interval_on_ref.0
                    + map_params.k as GnPosition
                    + 2 * c;
                _total_bases_contained_ref += int.interval_on_query.1 - int.interval_on_query.0
                    + map_params.k as GnPosition
                    + 2 * c;
            }

            let start =
                (i32::max(int.interval_on_query.0 as i32 - c as i32, 0)) as u32;
            let stop = int.interval_on_query.1 + c;
            all_intervals = all_intervals.union(&vec![(start, stop)].to_interval_set());
            //interval_vec.insert(int_insert, i);
            if sensitive_af{
                total_query_bases +=  int.query_range_len() - int.overlap + 2 * c + k as GnPosition;
                total_ref_range +=  int.query_range_len() - int.overlap + 2 * c + k as GnPosition;
            }

            avg_chain_int_len += int.query_range_len() - int.overlap + 2 * c + k as GnPosition;
            num_chains += 1;
        }

        if total_anchors == 0{
            continue;
        }

        if  total_range_query.1 - total_range_query.0 < map_params.min_length_cover as GnPosition{
            continue;
        }

        if !sensitive_af{
            total_query_bases += total_range_query.1 - total_range_query.0 + 2 * c + map_params.k as GnPosition;
            total_ref_range += total_range_query.1 - total_range_query.0 + 2 * c + map_params.k as GnPosition;
        }

        let mut num_seeds_in_intervals = 0;
        let mut upper_lower_seeds = 0;
        for pos in anchor_chunks.seeds_in_chunk[i].iter() {
            if all_intervals.contains(pos) {
                num_seeds_in_intervals += 1;
            }
        }

        let mut left_spacing_est =  0;
        let mut right_spacing_est = 0;
        let switched_ref_sketch;
        let switched_query_sketch;
        if switched{
            switched_query_sketch = &ref_sketch;
            switched_ref_sketch = &query_sketch;
        }
        else{
            switched_ref_sketch = &ref_sketch;
            switched_query_sketch = &query_sketch;
        }

        let ref_ctg_len = switched_ref_sketch.contig_lengths[leftmost_interval.ref_contig];
        trace!("switched ref ctg len {},id {}", ref_ctg_len, leftmost_interval.ref_contig);
        let q_ctg_len = switched_query_sketch.contig_lengths[leftmost_interval.query_contig];
        trace!("switched query ctg len {},id {}", q_ctg_len, leftmost_interval.query_contig);


        //TODO this was for testing... don't use this mechanism
        let extend = 0;
        if leftmost_interval.reverse_chain{
            let ref_ctg_len =switched_ref_sketch.contig_lengths[leftmost_interval.ref_contig];
            if ref_ctg_len - leftmost_interval.interval_on_ref.1 < extend{
                left_spacing_est = ref_ctg_len - leftmost_interval.interval_on_ref.1;
            }
        }
        else{
            //Clipped
            if leftmost_interval.interval_on_ref.0 < extend{
                left_spacing_est = leftmost_interval.interval_on_ref.0;
            }
        }
        if rightmost_interval.reverse_chain{
            if rightmost_interval.interval_on_ref.0 < extend{
                right_spacing_est = rightmost_interval.interval_on_ref.0;
            }
        }
        else{
            let ref_ctg_len = switched_ref_sketch.contig_lengths[rightmost_interval.ref_contig];
            trace!("{},{}",ref_ctg_len, rightmost_interval.ref_contig);
//                dbg!(ref_ctg_len, rightmost_interval);
//                dbg!(query_sketch.contig_lengths[leftmost_interval.ref_contig], leftmost_interval);
            if ref_ctg_len - rightmost_interval.interval_on_ref.1 < extend{
                right_spacing_est = ref_ctg_len - rightmost_interval.interval_on_ref.1;
            }
        }
//        dbg!(right_spacing_est, left_spacing_est);
        
        //        let spacing_est = (total_range_query.1 - total_range_query.0
        //            + 2 * ref_sketch.c as GnPosition)
        //            / num_seeds_in_intervals as GnPosition;
        for pos in anchor_chunks.seeds_in_chunk[i].iter() {
            if *pos + left_spacing_est >= total_range_query.0
                && *pos <= total_range_query.1 + right_spacing_est
            {
                upper_lower_seeds += 1;
            }
        }

        let mut anchors_in_chunk_considered = anchor_chunks.seeds_in_chunk[i].len();
        let putative_ani = f64::powf(
//                        total_anchors as f64 / (upper_lower_seeds) as f64,
            total_anchors as f64 / (num_seeds_in_intervals) as f64,
            1. / k as f64,
        );
        if putative_ani > 0.950
//            && total_bases_contained_query > ref_sketch.c as GnPosition * 20
            //&& total_bases_contained_query > c  * 3 * (upper_lower_seeds / total_anchors) as GnPosition
            && total_bases_contained_query > c * 4
            && !map_params.amino_acid
            && total_range_query.1 - total_range_query.0 < (CHUNK_SIZE_DNA * 9 / 10) as GnPosition 
            && anchors_in_chunk_considered as f64 > 1.05 * upper_lower_seeds as f64 
        {
            //                        anchors_in_chunk_considered = num_seeds_in_intervals;
            trace!("putative ani filter {} -> {}", anchors_in_chunk_considered, upper_lower_seeds);
            anchors_in_chunk_considered = upper_lower_seeds;
        }

        let test = false;
        if test{
            if right_spacing_est != 0 || left_spacing_est != 0{
                anchors_in_chunk_considered = upper_lower_seeds;
            }
            else{
                anchors_in_chunk_considered = anchor_chunks.seeds_in_chunk[i].len();
            }
        }

        let ml_hits = if map_params.amino_acid {
            f64::min(
                1.,
                total_anchors as f64 / anchors_in_chunk_considered as f64 * 6.,
            )
        } else {
            f64::min(
                1.,
                total_anchors as f64 / anchors_in_chunk_considered as f64,
            )
        };
        let ani_est = if map_params.amino_acid {
            f64::powf(ml_hits, 1. / k as f64)
        } else {
            f64::powf(ml_hits, 1. / k as f64)
        };

        //        total_bases_contained_query =
        //            total_range_query.1 - total_range_query.0 + map_params.k as GnPosition;
        //        total_bases_contained_ref =
        //            total_range_query.1 - total_range_query.0 + map_params.k as GnPosition;
        //
        //        total_query_range += total_bases_contained_query;
        //        total_ref_range += total_bases_contained_ref;

        //        ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i]));
        if map_params.amino_acid {
//            ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i].len() / 6));
            ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i].len() / 6));
//            ani_ests.push((ani_est, 1));
        } else {
            //ani_ests.push((ani_est, anchor_chunks.seeds_in_chunk[i].len()));
            ani_ests.push((ani_est, anchors_in_chunk_considered));
        }
        //                        ani_ests.push((ani_est, upper_lower_seeds));
        trace!(
            "Ani est fragment {}, total range {:?}, total anchors {}, seeds in fragment {:?},",
            ani_est,
            total_range_query,
            total_anchors,
            anchor_chunks.seeds_in_chunk[i].len(),
        );
        trace!(
            "Intervals {:?}, Num Anchors in Interval {}, Total Anchors {}",
            &intervals,
            intervals[0].num_anchors,
            total_anchors
        );
        _all_anchors_total += total_anchors;
        _num_good_chunks += 1;
    }
    ani_ests.sort_by(|x, y| x.partial_cmp(y).unwrap());

    if ani_ests.is_empty() || num_chains == 0 {
        let mut ret = AniEstResult::default();
        ret.ani = f32::NAN;
        return ret;
    }
    avg_chain_int_len /= num_chains;
    let mut total_multiplicitiy = 0;
    for ani in ani_ests.iter(){
        total_multiplicitiy += ani.1;
    }
    let lower;
    let upper;
    if map_params.median{
        lower = 0.499;
        upper = 0.501;
    } else if map_params.robust{
        lower = 0.10;
        upper = 0.90;
    } else {
        lower = 0.;
        upper = 1.;
    }

//    for ani_est in ani_ests.iter(){
//        println!("{},{}", ani_est.0, ani_est.1);
//    }

    let mut lower_i = 0;
    let mut upper_i = ani_ests.len()-1;
    let mut changed_l = false;
    let mut _changed_u = false;

    let mut curr_sum = 0;
    for (i,ani) in ani_ests.iter().enumerate(){
        curr_sum += ani.1;
        if curr_sum >= ((total_multiplicitiy as f64) * lower) as usize && !changed_l{
            lower_i = i;
            changed_l = true;
        }
        if curr_sum >= ((total_multiplicitiy as f64) * upper) as usize && !_changed_u{
            upper_i = i+1;
            _changed_u = true;
            break;
        }
    }

    let mut total_multiplicitiy = 0;
    let mut weighted_avg = 0.;
    for i in lower_i..upper_i{
        weighted_avg += ani_ests[i].0 * ani_ests[i].1 as f64;
        total_multiplicitiy += ani_ests[i].1;
        //        weighted_avg += ani_ests[i].0 * (ani_ests[i].1 as f64);
        //        total_weight_interval += ani_ests[i].1;
    }
    //    let mut final_ani = weighted_avg / total_weight_interval as f64;
    let mut final_ani = weighted_avg / total_multiplicitiy as f64;

//    let (upper, lower) = z_interval(&ani_ests);
    let ci_std = bootstrap_interval(&ani_ests);
    let ci = (ci_std.0, ci_std.1);
    let std = ci_std.2;
    let covered_query = f64::min(
        1.,
        total_query_bases as f64 / query_sketch.total_sequence_length as f64,
    );
    let covered_ref = f64::min(
        1.,
        total_ref_range as f64 / ref_sketch.total_sequence_length as f64,
    );
    
    let q_string = &query_sketch.file_name;
    let id_string = if map_params.amino_acid { "AAI" } else { "ANI" };
    trace!("Total 1-to-1 aligned bases: {}", total_query_bases);
    //println!("{}",total_query_bases);
    debug!(
        "Query {} Ref {} - {} {}, +/- = {}/{}. ",
        q_string,
        ref_sketch.file_name,
        id_string,
        final_ani,
        ci.0,
        ci.1,
    );

    if map_params.amino_acid{
        if covered_query < map_params.frac_cover_cutoff  || covered_ref < map_params.frac_cover_cutoff
        {
            final_ani = -1.;
        }
    }
    else if covered_query < map_params.frac_cover_cutoff  && covered_ref < map_params.frac_cover_cutoff
    {
        final_ani = -1.;
    }

    let mut sorted_contigs_q = query_sketch.contig_lengths.clone();
    let mut sorted_contigs_r = ref_sketch.contig_lengths.clone();
    sorted_contigs_q.sort();
    sorted_contigs_r.sort();
    let qs_lens = sorted_contigs_q.len();
    let rs_lens = sorted_contigs_r.len();
    let contig_quants_q = [sorted_contigs_q[qs_lens * 10 / 100 ], sorted_contigs_q[qs_lens* 50 / 100], sorted_contigs_q[qs_lens * 90 / 100]];
    let contig_quants_r = [sorted_contigs_r[rs_lens * 10 / 100 ], sorted_contigs_r[rs_lens* 50 / 100], sorted_contigs_r[rs_lens * 90 / 100]];
//    let mean_contig_len_q = query_sketch.contig_lengths.iter().map(|x| *x as f64).sum::<f64>()
//        /(query_sketch.contig_lengths.len() as f64);
//    let mean_contig_len_r = ref_sketch.contig_lengths.iter().map(|x| *x as f64).sum::<f64>()
//        /(ref_sketch.contig_lengths.len() as f64);

    AniEstResult {
        ani: final_ani as f32,
        align_fraction_query: covered_query as f32,
        align_fraction_ref: covered_ref as f32,
        ref_file: ref_sketch.file_name.clone(),
        query_file: query_sketch.file_name.clone(),
        query_contig: query_sketch.contigs[0].clone(),
        ref_contig: ref_sketch.contigs[0].clone(),
        num_contigs_r: ref_sketch.contigs.len() as u32,
        num_contigs_q: query_sketch.contigs.len() as u32,
        ci_upper: ci.1 as f32,
        ci_lower: ci.0 as f32,
        aai: map_params.amino_acid,
        quant_90_contig_len_q: contig_quants_q[2] as f32,
        quant_90_contig_len_r: contig_quants_r[2] as f32,
        quant_50_contig_len_q: contig_quants_q[1] as f32,
        quant_50_contig_len_r: contig_quants_r[1] as f32,
        quant_10_contig_len_q: contig_quants_q[0] as f32,
        quant_10_contig_len_r: contig_quants_r[0] as f32,
        std: std as f32,
        avg_chain_int_len,
        total_bases_covered: total_query_bases,
    }
}

#[inline]
pub fn score_anchors(anchor_curr: &Anchor, anchor_past: &Anchor, map_params: &MapParams) -> f64 {
    if anchor_curr.query_phase != anchor_past.query_phase
        || anchor_curr.ref_phase != anchor_past.ref_phase
    {
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

    if d_q > D_MAX_LIN_LENGTH || d_r > D_MAX_LIN_LENGTH {
        return f64::MIN;
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
    map_params.anchor_score - gap
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
    if ref_sketch.contig_lengths.is_empty() || query_sketch.contig_lengths.is_empty(){
        return (AnchorChunks::default(), true);
    }
//    let score_query = query_sketch.total_sequence_length as f64
//    let score_ref = ref_sketch.total_sequence_length as f64
//        * (ref_sketch.total_sequence_length as f64 / ref_sketch.contigs.len() as f64).ln();
//
    let mut ctgs_q = query_sketch.contig_lengths.iter().collect::<Vec<&GnPosition>>();
    let mean_ctg_len_q = query_sketch.contig_lengths.iter().map(|x| *x as f64).sum::<f64>()
        /(query_sketch.contig_lengths.len() as f64);
    ctgs_q.sort_unstable();
//    let med_ctg_len_q = *ctgs_q[query_sketch.contig_lengths.len()/2] as f64;
    let mean_ctg_len_r = ref_sketch.contig_lengths.iter().map(|x| *x as f64).sum::<f64>()
        /(ref_sketch.contig_lengths.len() as f64);

//    let score_query = (query_sketch.total_sequence_length as f64)
//        * f64::min(med_ctg_len_q, 40000.);
//    let score_ref = (ref_sketch.total_sequence_length as f64)
//        * f64::min(med_ctg_len_r, 40000.);

    let query_length_markers_proxy;
    let ref_length_markers_proxy;

    if query_sketch.total_sequence_length > 100_000 && ref_sketch.total_sequence_length > 100_000{
        query_length_markers_proxy = query_sketch.marker_seeds.len() as f64 * query_sketch.c as f64;
        ref_length_markers_proxy = ref_sketch.marker_seeds.len() as f64 * ref_sketch.c as f64;
    }
    else{
        query_length_markers_proxy = query_sketch.total_sequence_length as f64;
        ref_length_markers_proxy = ref_sketch.total_sequence_length as f64;
    }
    if switch_qr(mean_ctg_len_r,mean_ctg_len_q, query_length_markers_proxy, ref_length_markers_proxy, &query_sketch.file_name, &ref_sketch.file_name){
        switched = true;

        kmer_seeds_ref = query_sketch.kmer_seeds_k.as_ref().unwrap();
        kmer_seeds_query = ref_sketch.kmer_seeds_k.as_ref().unwrap();
        query_positions_all = vec![vec![]; ref_sketch.contigs.len()];
    } else {
        switched = false;

        kmer_seeds_ref = ref_sketch.kmer_seeds_k.as_ref().unwrap();
        kmer_seeds_query = query_sketch.kmer_seeds_k.as_ref().unwrap();
        query_positions_all = vec![vec![]; query_sketch.contigs.len()];
    }
    //    let kmer_seeds_ref = &ref_sketch.kmer_seeds_k[k];
    //    let kmer_seeds_query = &query_sketch.kmer_seeds_k[k];
    let mut anchors = vec![];
    let mut query_kmers_with_hits = 0;
    for (canon_kmer, query_pos) in kmer_seeds_query.iter() {
        if query_pos.len() > map_params.index_chain_band{
            continue;
        }
        let contains = kmer_seeds_ref.contains_key(canon_kmer);

        if !contains {
            for qpos in query_pos.iter() {
                query_positions_all[qpos.contig_index as usize].push(qpos.pos);
            }
        } else {
            let ref_pos = &kmer_seeds_ref[canon_kmer];

            if ref_pos.len() > map_params.index_chain_band{
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
        debug!(
            "no anchors found for {}, {}",
            &ref_sketch.file_name, &query_sketch.file_name
        );
        return (AnchorChunks::default(), true);
    }
    anchors.sort_unstable();
    for query_position_vec in query_positions_all.iter_mut() {
        query_position_vec.sort_unstable();
    }
    debug!(
        "Ref seeds len {}, Query seeds len {}, Anchors {}, Seeds hit query {}, Est {}, Ref_file {}, Query_file {}",
        kmer_seeds_ref.len(),
        kmer_seeds_query.len(),
        anchors.len(),
        query_kmers_with_hits,
        f64::powf(
            (query_kmers_with_hits as f64) / (kmer_seeds_query.len() as f64),
            1. / (k as f64)
        ),
        ref_sketch.file_name,
        query_sketch.file_name,
    );
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
            let mut _num_seeds_in_block = 0;
            let mut seed_pos_in_block = vec![];
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
                        trace!(
                            "start {}",
                            query_positions_all[last_query_contig as usize][running_counter]
                        );
                    }
                    seed_pos_in_block
                        .push(query_positions_all[last_query_contig as usize][running_counter]);
                    running_counter += 1;
                    _num_seeds_in_block += 1;
                } else {
                    trace!(
                        "end {}",
                        query_positions_all[last_query_contig as usize][running_counter - 1]
                    );
                    break;
                }
            }
            block_seeds.push(seed_pos_in_block);
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
    if !curr_anchor_chunk.is_empty() {
        let mut _num_seeds_in_block = 0;
        let mut seed_pos_in_block = vec![];
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
                seed_pos_in_block
                    .push(query_positions_all[last_query_contig as usize][running_counter]);
                running_counter += 1;
                _num_seeds_in_block += 1;
            } else {
                break;
            }
        }
        lengths.push(
            curr_anchor_chunk[curr_anchor_chunk.len() - 1].query_pos
                - curr_anchor_chunk[0].query_pos,
        );
        chunks.push(mem::take(&mut curr_anchor_chunk));
        block_seeds.push(seed_pos_in_block);
    }

    assert!(block_seeds.len() == chunks.len());

    (
        AnchorChunks {
            chunks,
            lengths,
            seeds_in_chunk: block_seeds,
        },
        switched,
    )
}

fn chain_anchors_ani(anchor_chunks: &AnchorChunks, map_params: &MapParams) -> Vec<ChainingResult> {
    let mut chaining_results = vec![];

    let num_chunks = anchor_chunks.chunks.len();
    let past_chain_length = usize::min(map_params.fragment_length / 2, map_params.bp_chain_band) as u32;

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
                if anchor_curr.ref_contig != anchor_past.ref_contig {
                    continue;
                }
                if anchor_curr.query_pos - anchor_past.query_pos > past_chain_length
                    || i - j > map_params.index_chain_band
                {
                    break;
                }
//                if anchor_curr.query_contig != anchor_past.query_contig {
//                    continue;
//                }
                if anchor_curr.ref_contig != anchor_past.ref_contig {
                    continue;
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
    chaining_results
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
        assert!(anchors[smallest_id].reverse_match == anchors[largest_id].reverse_match);
        assert!(anchors[smallest_id].ref_contig == anchors[largest_id].ref_contig);
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
            reverse_chain: anchors[smallest_id].reverse_match,
            overlap : 0
        };
        good_intervals.push(chain_interval);
    }
}
fn get_nonoverlapping_chains(
    intervals: &mut Vec<ChainInterval>,
    num_chunks: usize,
) -> Vec<Vec<ChainInterval>> {
    intervals.sort_by(|x, y| y.partial_cmp(x).unwrap());
    let mut interval_trees = FxHashMap::default();
    let mut interval_trees_ref = FxHashMap::default();
    let mut good_non_overlap_intervals = vec![vec![]; num_chunks];
    let mut bases_added = 0;
    for (i, int) in intervals.iter().enumerate() {
        let q_interval = (int.interval_on_query.0)..(int.interval_on_query.1);
        let r_interval = (int.interval_on_ref.0)..(int.interval_on_ref.1);
        let interval_tree_r = interval_trees_ref
            .entry(int.ref_contig)
            .or_insert(IntervalTree::new());
        let interval_tree_q = interval_trees
            .entry(int.query_contig)
            .or_insert(IntervalTree::new());

        let mut sum_overlaps_ref = 0;
        let mut sum_overlaps_query = 0;
        let no_overlap_ref;
        if interval_tree_r.find(&r_interval).count() == 0 {
            no_overlap_ref = true
        } else {
//            let mut TODO_intervals= vec![];
            let overlaps = interval_tree_r.find(&r_interval);
            let mut small_ol = false;
            for ol in overlaps {
                let ol_interval: &ChainInterval = &intervals[*ol.data()];
                let overlap = GnPosition::min(
                    int.interval_on_ref.1 - ol_interval.interval_on_ref.0,
                    ol_interval.interval_on_ref.1 - int.interval_on_ref.0,
                );
                sum_overlaps_ref += overlap;
//                TODO_intervals.push(ol_interval.clone());
                
            }
            if (sum_overlaps_ref as f32) < int.ref_range_len() as f32 * OVERLAP_ORTHOLOGOUS_FRACTION {
                bases_added += int.query_range_len();
                small_ol = true;
//                dbg!("ref", TODO_intervals, int);
            }
            if small_ol {
                no_overlap_ref = true;
            } else {
                no_overlap_ref = false;
            }
        }

        let no_overlap_query;
        if interval_tree_q.find(&q_interval).count() == 0 {
            no_overlap_query = true;
        } else {
            let overlaps = interval_tree_q.find(&q_interval);
            let mut small_ol = false;
//            let mut TODO_intervals= vec![];
            for ol in overlaps {
                let ol_interval: &ChainInterval = &intervals[*ol.data()];
                let overlap = GnPosition::min(
                    int.interval_on_query.1 - ol_interval.interval_on_query.0,
                    ol_interval.interval_on_query.1 - int.interval_on_query.0,
                );
                sum_overlaps_query += overlap;
//                TODO_intervals.push(ol_interval.clone());
            }

            if (sum_overlaps_query as f32) < int.query_range_len() as f32 * OVERLAP_ORTHOLOGOUS_FRACTION {
                bases_added += int.query_range_len();
                small_ol = true;
//                dbg!("query", TODO_intervals, int);
            }
            if small_ol {
                no_overlap_query = true;
            } else {
                no_overlap_query = false;
            }
        }
        if no_overlap_ref && no_overlap_query {
            //

            interval_tree_q.insert(q_interval, i);
            interval_tree_r.insert(r_interval, i);
            let mut cloned_int = int.clone();
            cloned_int.overlap = sum_overlaps_query;
            good_non_overlap_intervals[int.chunk_id].push(int.clone());
        }
    }
    trace!("Bases rescued by small overlapping orthologous threshold: {}", bases_added);
    //good_non_overlap_intervals.sort_by(|x, y| y.partial_cmp(&x).unwrap());
    good_non_overlap_intervals
}







































