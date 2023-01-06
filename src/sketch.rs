use crate::file_io;
use crate::params::*;
use crate::types::*;
use log::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;

use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

pub fn sketch(command_params: CommandParams, sketch_params: SketchParams) {
    let now = Instant::now();
    info!("Sketching files...");
    let p = format!("{}", command_params.out_file_name);
    if Path::new(&p).exists() {
        error!("Output directory exists; output directory must not be an existing directory. Exiting.");
        std::process::exit(1);
    }
    std::fs::create_dir_all(p).unwrap();

    let num_iters = command_params.ref_files.len();
    let mut shuffle_ref_files: Vec<&String> = command_params.ref_files.iter().collect();
    shuffle_ref_files.shuffle(&mut thread_rng());
    let counter: Mutex<usize> = Mutex::new(0);
    let marker_sketches: Mutex<Vec<Sketch>> = Mutex::new(vec![]);
    (0..num_iters).into_par_iter().for_each(|i| {
        let ref_sketches = file_io::fastx_to_sketches(
            //&vec![command_params.ref_files[i].clone()],
            &vec![shuffle_ref_files[i].clone()],
            &sketch_params,
            !command_params.screen,
        );
        let mut marker_ref_sketches = ref_sketches
            .iter()
            .map(|x| Sketch::get_markers_only(x))
            .collect::<Vec<Sketch>>();
        if ref_sketches.len() > 0 {
            let sketch = &ref_sketches[0];
            let marker_sketch = &mut marker_ref_sketches[0];
            let path = Path::new(&sketch.file_name);
            let filename = path.file_name().unwrap().to_str().unwrap();
            let mut file_bin = BufWriter::new(
                File::create(format!(
                    "{}/{}.sketch",
                    &command_params.out_file_name, filename
                ))
                .unwrap(),
            );

            trace!("{} compress factor", sketch.total_sequence_length / sketch.kmer_seeds_k.as_ref().unwrap().len());
            trace!("{} marker compress factor", sketch.total_sequence_length / sketch.marker_seeds.len());

            bincode::serialize_into(&mut file_bin, &(&sketch_params, sketch)).unwrap();

            let mut locked = marker_sketches.lock().unwrap();
            locked.push(std::mem::take(marker_sketch));
        }
        let mut locked = counter.lock().unwrap();
        *locked += 1;
        if *locked % 100 == 0 && *locked != 0 {
            info!("{} sequences sketched.", locked);
        }
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
    return;
}
