use crate::file_io;
use crate::params::*;
use crate::sketch_db::{SketchDbWriter};
use crate::types::*;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::sync::mpsc;

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
    std::fs::create_dir_all(&p).unwrap();

    if command_params.separate_sketches {
        if command_params.individual_contig_r {
            warn!("WARNING: --separate-sketches combined with -i (individual contigs) is NOT compatible with `skani search`. Use the default consolidated database format for search functionality with individual contigs.");
        }
        sketch_separate_files(command_params, sketch_params);
    } else {
        sketch_consolidated_db(command_params, sketch_params);
    }
    
    info!("Sketching time: {}", now.elapsed().as_secs_f32());
}

/// Create separate .sketch files (legacy format)
fn sketch_separate_files(command_params: CommandParams, sketch_params: SketchParams) {
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
        } else {
            ref_sketches = file_io::fastx_to_sketches(
                &vec![command_params.ref_files[i].clone()],
                &sketch_params,
                true,
            );
        }
        let marker_ref_sketches = ref_sketches
            .iter()
            .map(Sketch::get_markers_only)
            .collect::<Vec<Sketch>>();
            
        (0..ref_sketches.len()).into_par_iter().for_each(|j| {
            let sketch = &ref_sketches[j];
            let marker_sketch = &marker_ref_sketches[j];
            let path = Path::new(&sketch.file_name);
            let filename = path.file_name().unwrap().to_str().unwrap();
            let sketch_name;
            if command_params.individual_contig_r {
                sketch_name = format!(
                    "{}/{}_{}.sketch",
                    &command_params.out_file_name, j, filename
                );
            } else {
                sketch_name = format!(
                    "{}/{}.sketch",
                    &command_params.out_file_name, filename
                );
            }
            let mut file_bin = BufWriter::new(File::create(sketch_name).unwrap());

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
        File::create(format!("{}/markers.bin", &command_params.out_file_name)).unwrap(),
    );
    let markers = marker_sketches.into_inner().unwrap();
    bincode::serialize_into(&mut file_bin_marker, &(&sketch_params, markers)).unwrap();
}

/// Create consolidated database format using multi-producer single-consumer pattern
fn sketch_consolidated_db(command_params: CommandParams, sketch_params: SketchParams) {
    const CHANNEL_BUFFER_SIZE: usize = 1000; // Bounded channel to prevent memory blowup
    
    let (sender, receiver) = mpsc::sync_channel::<(Sketch, Sketch)>(CHANNEL_BUFFER_SIZE);
    let num_iters = command_params.ref_files.len();
    
    // Spawn consumer thread to handle sequential database writes
    let output_dir = command_params.out_file_name.clone();
    let sketch_params_consumer = sketch_params.clone();
    let consumer_handle = std::thread::spawn(move || {
        let mut db_writer = SketchDbWriter::new(&output_dir)
            .unwrap_or_else(|e| {
                error!("Failed to create consolidated database writer: {}", e);
                std::process::exit(1);
            });
        
        let mut marker_sketches = Vec::new();
        let mut sketch_count = 0;
        
        // Consume sketches from the channel
        for (sketch, marker_sketch) in receiver {
            trace!("{} compress factor", sketch.total_sequence_length / sketch.kmer_seeds_k.as_ref().unwrap().len());
            trace!("{} marker compress factor", sketch.total_sequence_length / sketch.marker_seeds.len());

            // Add sketch to consolidated database
            if let Err(e) = db_writer.add_sketch(&sketch_params_consumer, &sketch) {
                error!("Failed to add sketch to database: {}", e);
                std::process::exit(1);
            }

            marker_sketches.push(marker_sketch);
            sketch_count += 1;
            
            if sketch_count % 100 == 0 && sketch_count != 0 {
                info!("{} sequences sketched.", sketch_count);
            }
        }

        // Finalize the consolidated database
        if let Err(e) = db_writer.finalize(&output_dir) {
            error!("Failed to finalize consolidated database: {}", e);
            std::process::exit(1);
        }

        // Write markers.bin (still separate for compatibility)
        let mut file_bin_marker = BufWriter::new(
            File::create(format!("{}/markers.bin", &output_dir)).unwrap(),
        );
        bincode::serialize_into(&mut file_bin_marker, &(&sketch_params_consumer, marker_sketches)).unwrap();
        
        sketch_count
    });

    // Multi-producer: parallel sketch generation
    (0..num_iters).into_par_iter().for_each(|i| {
        let sender = sender.clone();
        let ref_sketches;
        
        if command_params.individual_contig_r {
            ref_sketches = file_io::fastx_to_multiple_sketch_rewrite(
                &vec![command_params.ref_files[i].clone()],
                &sketch_params,
                true,
            );
        } else {
            ref_sketches = file_io::fastx_to_sketches(
                &vec![command_params.ref_files[i].clone()],
                &sketch_params,
                true,
            );
        }

        for sketch in ref_sketches {
            let marker_sketch = Sketch::get_markers_only(&sketch);
            
            // Send to consumer thread (this will block if channel is full)
            if sender.send((sketch, marker_sketch)).is_err() {
                error!("Failed to send sketch to consumer thread");
                std::process::exit(1);
            }
        }
    });

    // Drop sender to signal end of production
    drop(sender);
    
    // Wait for consumer to finish
    match consumer_handle.join() {
        Ok(final_count) => {
            info!("Successfully wrote {} sketches to consolidated database", final_count);
        }
        Err(_) => {
            error!("Consumer thread panicked");
            std::process::exit(1);
        }
    }
}
