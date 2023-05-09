use crate::file_io;
use crate::params::*;
use crate::types::*;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;

use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

pub fn sketch(command_params: CommandParams, sketch_params: SketchParams) {
    let now = Instant::now();
    info!("Sketching files...");
    let p = command_params.out_file_name.to_string();
    if Path::new(&p).exists() {
        error!("Output directory exists; output directory must not be an existing directory. Exiting.");
        std::process::exit(1);
    }
    std::fs::create_dir_all(p).unwrap();

    let num_iters = command_params.ref_files.len();
    let counter: Mutex<usize> = Mutex::new(0);
    let marker_sketches: Mutex<Vec<Sketch>> = Mutex::new(vec![]);
    (0..num_iters).into_par_iter().for_each(|i| {
        let ref_sketches;
        if command_params.individual_contig_r {
            ref_sketches = file_io::fastx_to_multiple_sketch_rewrite(
                &vec![command_params.ref_files[i].clone()],
                &sketch_params,
                true,
        );
        }
        else{
            ref_sketches = file_io::fastx_to_sketches(
                &vec![command_params.ref_files[i].clone()],
                &sketch_params,
                true,);
        }
        let marker_ref_sketches = ref_sketches
            .iter()
            .map(Sketch::get_markers_only)
            .collect::<Vec<Sketch>>();
        (0..ref_sketches.len()).into_par_iter().for_each(|j|{
            let sketch = &ref_sketches[j];
            let marker_sketch = &marker_ref_sketches[j];
            let path = Path::new(&sketch.file_name);
            let filename = path.file_name().unwrap().to_str().unwrap();
            let sketch_name;
            if command_params.individual_contig_r{
                sketch_name = format!(
                    "{}/{}_{}.sketch",
                    &command_params.out_file_name, j, filename
                );
            }
            else{
                sketch_name = format!(
                    "{}/{}.sketch",
                    &command_params.out_file_name, filename
                );
            }
            let mut file_bin = BufWriter::new(
                File::create(sketch_name)
                .unwrap(),
            );

            trace!("{} compress factor", sketch.total_sequence_length / sketch.kmer_seeds_k.as_ref().unwrap().len());
            trace!("{} marker compress factor", sketch.total_sequence_length / sketch.marker_seeds.len());

            bincode::serialize_into(&mut file_bin, &(&sketch_params, sketch)).unwrap();

            let mut locked = marker_sketches.lock().unwrap();
            locked.push(marker_sketch.clone());
            let mut locked = counter.lock().unwrap();
            *locked += 1;
            if *locked % 100 == 0 && *locked != 0 {
                info!("{} sequences sketched.", locked);
            }
        });
    });
    let mut file_bin_marker = BufWriter::new(
        File::create(format!(
            "{}/markers.bin",
            &command_params.out_file_name
        ))
        .unwrap(),
    );
    let markers = marker_sketches.into_inner().unwrap();
    bincode::serialize_into(&mut file_bin_marker, &(&sketch_params, markers)).unwrap();
    info!("Sketching time: {}", now.elapsed().as_secs_f32());
    
}
