use crate::file_io;
use crate::types::*;
use std::fs::File;
use std::path::Path;
use crate::params::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rayon::prelude::*;
use std::io::{BufWriter};

pub fn sketch(command_params: CommandParams, sketch_params: SketchParams) {
    let p = format!("{}", command_params.out_file_name);
    std::fs::create_dir_all(p).unwrap();

    let num_iters = command_params.ref_files.len();
    let mut shuffle_ref_files: Vec<&String> = command_params.ref_files.iter().collect();
    shuffle_ref_files.shuffle(&mut thread_rng());
    (0..num_iters).into_par_iter().for_each(|i| {
        let ref_sketches = file_io::fastx_to_sketches(
            //&vec![command_params.ref_files[i].clone()],
            &vec![shuffle_ref_files[i].clone()],
            &sketch_params,
            !command_params.screen,
        );
        let marker_ref_sketches = ref_sketches
            .iter()
            .map(|x| Sketch::get_markers_only(x))
            .collect::<Vec<Sketch>>();
        if ref_sketches.len() > 0 {
            let sketch = &ref_sketches[0];
            let marker_sketch = &marker_ref_sketches[0];
            let path = Path::new(&sketch.file_name);
            let filename = path.file_name().unwrap().to_str().unwrap();
            let mut file_bin = BufWriter::new(
                File::create(format!(
                    "{}/{}.sketch",
                    &command_params.out_file_name, filename
                ))
                .unwrap(),
            );
            let mut file_bin_marker = BufWriter::new(
                File::create(format!(
                    "{}/{}.marker",
                    &command_params.out_file_name, filename
                ))
                .unwrap(),
            );

            bincode::serialize_into(&mut file_bin, &(&sketch_params, sketch)).unwrap();

            bincode::serialize_into(&mut file_bin_marker, &(&sketch_params, marker_sketch))
                .unwrap();
        }
    });
    return;
}
