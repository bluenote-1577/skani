use crate::file_io;
use crate::params::*;
use crate::sketch_db::{SketchDbWriter, INDEX_DB, MARKERS_FILE, MARKERS_FILE_LEGACY, MARKERS_TOMBSTONE, SKETCHES_DB};
use crate::types::*;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::sync::mpsc;

use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

/// Drop a readable stub at the pre-0.4 `markers.bin` name. See
/// `sketch_db::MARKERS_TOMBSTONE` for why this exists.
fn write_markers_tombstone(output_dir: &str) {
    let path = format!("{}/{}", output_dir, MARKERS_FILE_LEGACY);
    if let Err(e) = std::fs::write(&path, MARKERS_TOMBSTONE) {
        warn!("Could not write {}: {}", path, e);
    }
}

/// Days since the Unix epoch -> (year, month, day). Howard Hinnant's
/// civil_from_days; avoids pulling in a date crate for one timestamp.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn utc_timestamp() -> String {
    match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        Ok(dur) => {
            let secs = dur.as_secs() as i64;
            let (y, m, d) = civil_from_days(secs.div_euclid(86400));
            let tod = secs.rem_euclid(86400);
            format!(
                "{:04}-{:02}-{:02} {:02}:{:02}:{:02} UTC",
                y, m, d, tod / 3600, (tod % 3600) / 60, tod % 60
            )
        }
        Err(_) => "unknown".to_string(),
    }
}

/// Write a README.txt next to the sketches recording how the database was
/// built and what to do when something does not work. Purely informational -
/// skani never reads it back.
fn write_database_readme(
    command_params: &CommandParams,
    sketch_params: &SketchParams,
    separate_sketches: bool,
    num_genomes: usize,
) {
    let dir = &command_params.out_file_name;
    let layout = if separate_sketches {
        format!(
            "  *.sketch       one sketch per input genome\n  {}   marker k-mers, used to screen candidates before full comparison\n  {}    compatibility stub, see below\n  README.txt     this file",
            MARKERS_FILE, MARKERS_FILE_LEGACY
        )
    } else {
        format!(
            "  {}   the sketches, one record per genome\n  {}      byte offsets into {}\n  {}   marker k-mers, used to screen candidates before full comparison\n  {}    compatibility stub, see below\n  README.txt     this file",
            SKETCHES_DB, INDEX_DB, SKETCHES_DB, MARKERS_FILE, MARKERS_FILE_LEGACY
        )
    };

    let text = format!(
"skani database
==============

Built by   : skani {version}
Built on   : {when}
Command    : {cmd}

Sketch parameters (baked in -- see the last note)
-------------------------------------------------
  c (compression)       : {c}
  k (k-mer length)      : {k}
  marker c              : {marker_c}
  sequence type         : {seqtype}
  individual contigs    : {indiv}
  inputs sketched       : {n}

Files
-----
{layout}

Notes
-----
* This database is in the skani 0.4 format and CANNOT be read by skani 0.3.x
  or earlier. Point an older skani at this folder and it exits with
  'Problem reading .../{markers_legacy}'. That file is a deliberate stub whose
  only job is to make old versions fail immediately rather than slowly; the
  real marker data is in {markers}.

* skani 0.4 still reads databases built by 0.3.x. It warns once and converts
  them while loading, which is slower than using a 0.4 database.

* Search against this database with:
      skani search -d {dir} -q query.fa

* The sketch parameters above are fixed at build time. Comparing against
  genomes sketched with a different -c or -k does not work -- if you see a
  warning about 'Input parameter c = ... is not equal to the sketch
  parameter c = ...', that is the cause, and the sketch parameters win.
",
        version = env!("CARGO_PKG_VERSION"),
        when = utc_timestamp(),
        cmd = std::env::args().collect::<Vec<String>>().join(" "),
        c = sketch_params.c,
        k = sketch_params.k,
        marker_c = sketch_params.marker_c,
        seqtype = if sketch_params.use_aa { "amino acid (--aai)" } else { "DNA" },
        indiv = if command_params.individual_contig_r { "yes (-i)" } else { "no" },
        n = num_genomes,
        layout = layout,
        markers = MARKERS_FILE,
        markers_legacy = MARKERS_FILE_LEGACY,
        dir = dir,
    );

    let path = format!("{}/README.txt", dir);
    if let Err(e) = std::fs::write(&path, text) {
        warn!("Could not write {}: {}", path, e);
    }
}

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

            trace!("{} compress factor", sketch.total_sequence_length / sketch.num_distinct_seeds());
            trace!("{} marker compress factor", sketch.total_sequence_length / sketch.marker_seeds.len());

            write_sketch_header(&mut file_bin).unwrap();
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
        File::create(format!("{}/{}", &command_params.out_file_name, MARKERS_FILE)).unwrap(),
    );
    let markers = marker_sketches.into_inner().unwrap();
    write_sketch_header(&mut file_bin_marker).unwrap();
    bincode::serialize_into(&mut file_bin_marker, &(&sketch_params, markers)).unwrap();
    write_markers_tombstone(&command_params.out_file_name);
    write_database_readme(&command_params, &sketch_params, true, num_iters);
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
            trace!("{} compress factor", sketch.total_sequence_length / sketch.num_distinct_seeds());
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

        // Write the marker file, plus a tombstone at the pre-0.4 name.
        let mut file_bin_marker = BufWriter::new(
            File::create(format!("{}/{}", &output_dir, MARKERS_FILE)).unwrap(),
        );
        write_sketch_header(&mut file_bin_marker).unwrap();
        bincode::serialize_into(&mut file_bin_marker, &(&sketch_params_consumer, marker_sketches)).unwrap();
        write_markers_tombstone(&output_dir);
        
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
            write_database_readme(&command_params, &sketch_params, false, final_count);
        }
        Err(_) => {
            error!("Consumer thread panicked");
            std::process::exit(1);
        }
    }
}
