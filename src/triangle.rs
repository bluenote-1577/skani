use crate::chain;
use crate::file_io;
use crate::params::*;
use crate::screen;
use crate::types::*;
use fxhash::{FxHashMap};
use log::*;
use rayon::prelude::*;
use std::sync::Mutex;
use std::time::Instant;

pub fn triangle(command_params: CommandParams, mut sketch_params: SketchParams) {
    let ref_sketches;
    let now = Instant::now();
    if command_params.refs_are_sketch {
        info!("Sketches detected.");
        (sketch_params, ref_sketches) = file_io::sketches_from_sketch(
            &command_params.ref_files,
        );
    } else if command_params.individual_contig_r {
        ref_sketches = file_io::fastx_to_multiple_sketch_rewrite(
            &command_params.ref_files,
            &sketch_params,
            true,
        );
    } else {
        ref_sketches =
            file_io::fastx_to_sketches(&command_params.ref_files, &sketch_params, true);
    }
    let screen_val;
    if command_params.screen_val == 0. {
        if sketch_params.use_aa {
            screen_val = SEARCH_AAI_CUTOFF_DEFAULT;
        } else {
            screen_val = SEARCH_ANI_CUTOFF_DEFAULT;
        }
    } else {
        screen_val = command_params.screen_val;
    }
    let anis: Mutex<FxHashMap<usize, FxHashMap<usize, AniEstResult>>> =
        Mutex::new(FxHashMap::default());

    if ref_sketches.is_empty() {
        error!("No genomes/sketches found.");
        std::process::exit(1)
    }
    let kmer_to_sketch = screen::kmer_to_sketch_from_refs(&ref_sketches);
    let counter : Mutex<usize> = Mutex::new(0);

    (0..ref_sketches.len() - 1)
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
            
            let ref_sketch_i = &ref_sketches[i];
            //if command_params.screen {
            let screened_refs = screen::screen_refs(
                screen_val,
                &kmer_to_sketch,
                ref_sketch_i,
                &sketch_params,
                &ref_sketches,
            );
            debug!(
                "{} has {} refs passing screening.",
                ref_sketch_i.file_name,
                screened_refs.len()
            );

            screened_refs.into_par_iter().for_each(|j| {
                if j > i {
                    let map_params = chain::map_params_from_sketch(
                        ref_sketch_i,
                        sketch_params.use_aa,
                        &command_params,
                    );
                    let ref_sketch_j = &ref_sketches[j];
                    let ani_res = chain::chain_seeds(ref_sketch_i, ref_sketch_j, map_params);
                    let mut locked = anis.lock().unwrap();
                    let mapi = locked.entry(i).or_insert(FxHashMap::default());
                    mapi.insert(j, ani_res);
                }
            });
            let mut locked = counter.lock().unwrap();
            *locked += 1;
            if *locked % 100 == 0 && *locked != 0{
                info!("{} query sequences processed.", locked);
            }
        });
    let anis = anis.into_inner().unwrap();
    if command_params.sparse {
        file_io::write_sparse_matrix(&anis, &ref_sketches, &command_params.out_file_name, sketch_params.use_aa, command_params.est_ci);
    } else {
        file_io::write_phyllip_matrix(
            &anis,
            &ref_sketches,
            &command_params.out_file_name,
            command_params.individual_contig_r,
            command_params.full_matrix,
            sketch_params.use_aa
        );
    }
    info!("ANI/AAI triangle time: {}", now.elapsed().as_secs_f32());
}
