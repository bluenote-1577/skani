use crate::types::*;
use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;

pub fn predict_from_ani_res(ani_res: &mut AniEstResult, model: &GBDT) {
    if ani_res.ani > 0.9 {
        let data;
        if ani_res.quant_50_contig_len_r > ani_res.quant_50_contig_len_q {
            data = Data::new_test_data(
                vec![
                    ani_res.ani * 100.,
                    ani_res.align_fraction_ref * 100.,
                    ani_res.align_fraction_query * 100.,
                    ani_res.ci_lower * 100.,
                    ani_res.ci_upper * 100.,
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
                    ani_res.align_fraction_ref * 100.,
                    ani_res.align_fraction_query * 100.,
                    ani_res.ci_lower * 100.,
                    ani_res.ci_upper * 100.,
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
            ani_res.ani = pred_ani_res / 100.;
        }
    }
}
