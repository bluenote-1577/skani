use crate::types::*;
use crate::model;
use crate::params::*;
use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;
use log::*;

pub fn get_model(c: usize, learned_ani: bool) -> Option<GBDT>{
    let model: Option<GBDT>;
    if learned_ani {
        if (c as i32 - 125 as i32).abs() < (c as i32 - 200 as i32).abs(){
            debug!("Using C125 regression model.");
            model = Some(serde_json::from_str(model::MODEL).unwrap());
        }
        else{
            debug!("Using C200 regression model.");
            model = Some(serde_json::from_str(model::MODEL_C200).unwrap());
        }
    }
    else{
        model = None;
    }
    model
}

pub fn predict_from_ani_res(ani_res: &mut AniEstResult, model: &GBDT) {
    if ani_res.ani > 0.9  && ani_res.total_bases_covered > TOTAL_BASES_REGRESS_CUTOFF as GnPosition{
        let data;
        if ani_res.quant_50_contig_len_r > ani_res.quant_50_contig_len_q {
            data = Data::new_test_data(
                vec![
                    ani_res.ani * 100.,
                    ani_res.align_fraction_ref * 100.,
                    ani_res.align_fraction_query * 100.,
                    ani_res.std,
                    ani_res.quant_90_contig_len_r,
                    ani_res.quant_50_contig_len_r,
                    ani_res.quant_10_contig_len_r,
                    ani_res.quant_90_contig_len_q,
                    ani_res.quant_50_contig_len_q,
                    ani_res.quant_10_contig_len_q,
                ],
                None,
            );
        } else {
            data = Data::new_test_data(
                vec![
                    ani_res.ani * 100.,
                    ani_res.align_fraction_query * 100.,
                    ani_res.align_fraction_ref * 100.,
                    ani_res.std,
                    ani_res.quant_90_contig_len_q,
                    ani_res.quant_50_contig_len_q,
                    ani_res.quant_10_contig_len_q,
                    ani_res.quant_90_contig_len_r,
                    ani_res.quant_50_contig_len_r,
                    ani_res.quant_10_contig_len_r,
                ],
                None,
            );
        }
        let pred_ani_res = model.predict(&vec![data])[0];
        //dbg!(ani_res.ani_res*100., pred_ani_res);
        if pred_ani_res < 100. {
            ani_res.ci_upper = (ani_res.ci_upper - ani_res.ani) + pred_ani_res / 100.;
            ani_res.ci_lower = (ani_res.ci_lower - ani_res.ani) + pred_ani_res / 100.;
            ani_res.ani = pred_ani_res / 100.;
        }
    }
}
