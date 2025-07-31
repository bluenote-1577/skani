use skani::chain::*;
use skani::seeding::*;
use skani::avx2_seeding::*;
use skani::regression::*;
use skani::file_io::*;
use skani::params::*;
use skani::types::*;
fn default_params(mode: Mode) -> (CommandParams, SketchParams) {
    let cmd_params = CommandParams {
        screen: false,
        screen_val: 0.00,
        mode: mode,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        diagonal: false,
        max_results: 10000000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: 0.15,
        keep_refs: false,
        est_ci: false,
        learned_ani: true,
        detailed_out: false,
        distance: false,
        rescue_small: true,
        separate_sketches: false
    };

    let sketch_params = SketchParams::new(1000, 125, 15, false, false);
    return (cmd_params, sketch_params);
}

#[test]
fn fast_ecoli_test_simple() {
    let (mut command_params, sketch_params) = default_params(Mode::Dist);
    command_params
        .ref_files
        .push("./test_files/e.coli-W.fasta.gz".to_string());
    command_params
        .query_files
        .push("./test_files/e.coli-W.fasta.gz".to_string());

    let ref_sketch = fastx_to_sketches(&command_params.ref_files, &sketch_params, true)[0].clone();
    let query_sketch =
        fastx_to_sketches(&command_params.query_files, &sketch_params, true)[0].clone();
    let map_params = map_params_from_sketch(&ref_sketch, sketch_params.use_aa, &command_params, &None);
    let ani_res = chain_seeds(&ref_sketch, &query_sketch, map_params);
    assert!(ani_res.ani >= 1.0);
    assert!(ani_res.align_fraction_query >= 0.99);
    assert!(ani_res.align_fraction_ref >= 0.99);
}

#[test]
fn fast_ecoli_plasmid_test() {
    let (mut command_params, sketch_params) = default_params(Mode::Dist);
    command_params
        .ref_files
        .push("./test_files/e.coli-o157.fasta.sketch".to_string());
    command_params
        .query_files
        .push("./test_files/o157_plasmid.fasta".to_string());

    let ref_sketch = sketches_from_sketch(&command_params.ref_files).1[0].clone();
    let query_sketch =
        fastx_to_sketches(&command_params.query_files, &sketch_params, true)[0].clone();
    let map_params = map_params_from_sketch(&ref_sketch, sketch_params.use_aa, &command_params, &None);
    let ani_res = chain_seeds(&ref_sketch, &query_sketch, map_params);
    assert!(ani_res.ani >= 1.0);
    assert!(ani_res.align_fraction_query >= 0.99);
    assert!(ani_res.align_fraction_ref >= 0.005);
}

#[test]
fn fast_eukaryote_test() {

/*   Mummer stats
    [Bases]
    TotalBases                  37074820             45684866
    AlignedBases      28626642(77.2132%)   28594130(62.5899%)
    UnalignedBases     8448178(22.7868%)   17090736(37.4101%)

    [Alignments]
    1-to-1                         10157                10157
    TotalLength                 28580162             28584588
    AvgLength                  2813.8389            2814.2747
    AvgIdentity                  98.6310              98.6310
*/

    let (mut command_params, sketch_params) = default_params(Mode::Dist);
    command_params
        .ref_files
        .push("./test_files/TOPAZ_IOD1_E001.fna.gz".to_string());
    command_params
        .query_files
        .push("./test_files/TOPAZ_RSS1_E007.fna.gz".to_string());

    let ref_sketch = fastx_to_sketches(&command_params.ref_files, &sketch_params, true)[0].clone();
    let query_sketch =
        fastx_to_sketches(&command_params.query_files, &sketch_params, true)[0].clone();
    let map_params = map_params_from_sketch(&ref_sketch, sketch_params.use_aa, &command_params, &None);
    let ani_res = chain_seeds(&ref_sketch, &query_sketch, map_params);

    assert!(ani_res.ani >= 0.98);
    assert!(ani_res.align_fraction_ref>= 0.70);
    assert!(ani_res.align_fraction_query >= 0.55);

    assert!(ani_res.align_fraction_ref<= 0.80);
    assert!(ani_res.align_fraction_query <= 0.65);
    let old_ani = ani_res.ani;


    let model_opt = get_model(sketch_params.c, command_params.learned_ani);
    let map_params = map_params_from_sketch(&ref_sketch, sketch_params.use_aa, &command_params, &model_opt);
    let ani_res = chain_seeds(&ref_sketch, &query_sketch, map_params);

    assert!(ani_res.ani >= 0.98);
    assert!(ani_res.ani <= old_ani);
    println!("{:?}", ani_res);
}

#[test]
fn fast_avx2_vs_normal_code(){
    let str1 = b"ATCAGATTTAAAAAAAAATTTTGCTAGCTGATCGATCGATCGATGTGTATATATTAAAAGAGAGAGAGGGGGGGGAAAAAAAAAAAAACTGATCGATCGATGCTAGCTAGTCAGTCGATG";
    let (_command_params, mut sketch_params) = default_params(Mode::Dist);
    sketch_params.c = 10;
    let mut new_sketch1 = Sketch::default();
    let mut new_sketch2 = Sketch::default();
    unsafe{
        avx2_fmh_seeds(str1, &sketch_params, 0, &mut new_sketch1, true);
    }

    fmh_seeds(str1, &sketch_params, 0, &mut new_sketch2, true);
    assert!(new_sketch1 == new_sketch2);
    //println!("{:?}", new_sketch1.kmer_seeds_k.unwrap());
}


//Ns are treated as As right now. Luckily, the hash function doesn't accept AAAAAA... 
//for FMH when c = 30. 
#[test]
fn fast_NNN_test_code(){
    let str1 = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNn";
    let (mut _command_params, mut sketch_params) = default_params(Mode::Dist);
    sketch_params.c = 30;
    let mut new_sketch1 = Sketch::default();
    fmh_seeds(str1, &sketch_params, 0, &mut new_sketch1, true);
    assert!(new_sketch1.kmer_seeds_k.unwrap().len() == 0);
}
