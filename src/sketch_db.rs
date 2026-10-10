use crate::params::*;
use crate::types::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use log::*;

/// Index entry for locating sketches in the concatenated database
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct IndexEntry {
    pub file_name: String,
    pub offset: u64,
    pub length: u64,
}

/// Writer for creating consolidated sketch databases
/// Filenames for the 0.4 consolidated database. The 0.4 record layout is not
/// readable by skani <= 0.3.x, so the files are named distinctly: an older
/// binary pointed at a 0.4 database simply does not find the names it expects,
/// instead of finding them and failing to parse them.
///
/// `markers.bin` deliberately keeps its name. An older skani that finds it and
/// fails to parse it reports "Problem reading markers.bin", which says more
/// than the "markers.bin not found" it would report if that file were renamed
/// too.
/// The 0.4 marker file, and the pre-0.4 name it replaces.
pub const DB_README: &str = "README.txt";
pub const MARKERS_FILE: &str = "markers.bin2";
pub const MARKERS_FILE_LEGACY: &str = "markers.bin";

/// Written at `markers.bin` in a 0.4 database so that skani <= 0.3.x fails
/// immediately instead of grinding through a format it cannot read.
///
/// `search` on an older binary reads `markers.bin` first and deserializes it
/// as `(SketchParams, Vec<Sketch>)`. When that file is a real 0.4 marker file,
/// the old parser consumes the header as `c`/`k`/`marker_c`, finds plausible
/// values for the two bools that follow, and keeps going -- so it does a lot of
/// work before failing. Byte 24 of the stream is `SketchParams::use_syncs`, a
/// bool, and bincode rejects any byte other than 0 or 1. Every byte of this
/// message is printable ASCII (>= 0x20), so the old parser is guaranteed to
/// stop after 25 bytes. It is also readable if someone simply cats the file.
pub const MARKERS_TOMBSTONE: &str = "\
SKANI-DB-V2: this database requires skani >= 0.4.0.
The marker data lives in markers.bin2. This file exists only so that older
versions of skani fail immediately with a clear error instead of slowly
parsing a format they cannot read.
Upgrade skani, or rebuild this database with skani >= 0.4.0.
";

pub const SKETCHES_DB: &str = "sketches.db2";
pub const INDEX_DB: &str = "index.db2";
/// Pre-0.4 names, still read.
pub const SKETCHES_DB_LEGACY: &str = "sketches.db";
pub const INDEX_DB_LEGACY: &str = "index.db";

pub struct SketchDbWriter {
    concat_file: BufWriter<File>,
    index: Vec<IndexEntry>,
    current_offset: u64,
}

/// Reader for accessing consolidated sketch databases
pub struct SketchDbReader {
    mmap: memmap2::Mmap,
    index: Vec<(u64, u64)>, // Vector of (offset, length) pairs, indexed by sketch index
}

impl SketchDbWriter {
    /// Create a new consolidated sketch database writer
    pub fn new(output_dir: &str) -> Result<Self, Box<dyn std::error::Error>> {
                
        let concat_path = format!("{}/{}", output_dir, SKETCHES_DB);
        let concat_file = BufWriter::new(File::create(concat_path)?);
        
        Ok(SketchDbWriter {
            concat_file,
            index: Vec::new(),
            current_offset: 0,
        })
    }

    /// Add a sketch to the consolidated database
    pub fn add_sketch(&mut self, sketch_params: &SketchParams, sketch: &Sketch) -> Result<(), Box<dyn std::error::Error>> {
        // Serialize the sketch data. Each record carries its own format header,
        // exactly like a standalone .sketch file, so a record can be classified
        // from its own bytes without consulting anything else in the database.
        // The header is part of the record, so index offsets and lengths cover it.
        let mut serialized = sketch_header_bytes().to_vec();
        bincode::serialize_into(&mut serialized, &(sketch_params, sketch))?;
        let length = serialized.len() as u64;

        // Record the index entry
        let entry = IndexEntry {
            file_name: sketch.file_name.clone(),
            offset: self.current_offset,
            length,
        };
        self.index.push(entry);

        // Write to the main database file
        self.concat_file.write_all(&serialized)?;
        self.current_offset += length;

        trace!("Added sketch {} at offset {} with length {}", sketch.file_name, self.current_offset - length, length);
        Ok(())
    }

    /// Finalize the database by writing the index
    pub fn finalize(mut self, output_dir: &str) -> Result<(), Box<dyn std::error::Error>> {
        // Flush and close main database file
        self.concat_file.flush()?;
        drop(self.concat_file);

        // Write index.db
        let index_path = format!("{}/{}", output_dir, INDEX_DB);
        let index_file = File::create(index_path)?;
        let mut index_writer = BufWriter::new(index_file);
        bincode::serialize_into(&mut index_writer, &self.index)?;
        index_writer.flush()?;

        info!("Consolidated sketch database written with {} sketches", self.index.len());
        Ok(())
    }
}

impl SketchDbReader {
    /// Open a consolidated sketch database for reading
    pub fn new(database_dir: &str) -> Result<Self, Box<dyn std::error::Error>> {
        // Load index.db
        // Prefer the 0.4 names; fall back to the pre-0.4 ones. A database
        // with the old names always holds old-format records, and those are
        // detected and warned about per record on deserialization.
        let v2_index = format!("{}/{}", database_dir, INDEX_DB);
        let (index_path, concat_name) = if Path::new(&v2_index).exists() {
            (v2_index, SKETCHES_DB)
        } else {
            (format!("{}/{}", database_dir, INDEX_DB_LEGACY), SKETCHES_DB_LEGACY)
        };
        let index_file = File::open(index_path)?;
        let index_reader = BufReader::new(index_file);
        let index_vec: Vec<IndexEntry> = bincode::deserialize_from(index_reader)?;

        // Convert to vector for index-based lookups
        let index: Vec<(u64, u64)> = index_vec.iter()
            .map(|entry| (entry.offset, entry.length))
            .collect();

        // Memory map the main database file
        let concat_path = format!("{}/{}", database_dir, concat_name);
        let concat_file = File::open(concat_path)?;
        let mmap = unsafe { memmap2::Mmap::map(&concat_file)? };

        info!("Loaded consolidated sketch database with {} sketches", index.len());
        Ok(SketchDbReader { mmap, index })
    }

    /// Get a sketch by index
    pub fn get_sketch(&self, index: usize) -> Result<(SketchParams, Sketch), Box<dyn std::error::Error>> {
        if let Some(&(offset, length)) = self.index.get(index) {
            let start = offset as usize;
            let end = start + length as usize;
            
            // Time the mmap slice access (potential disk I/O)
            let bytes = &self.mmap[start..end];
            
            // Time the deserialization (memory copy + parsing). Each record is
            // self-describing: a current one opens with the format header, a
            // pre-0.4 one opens straight into SketchParams.
            match detect_sketch_format(bytes) {
                SketchFormat::Native => {
                    let (params, sketch) = bincode::deserialize(&bytes[SKETCH_HEADER_LEN..])?;
                    Ok((params, sketch))
                }
                SketchFormat::Legacy => {
                    warn_legacy_sketch_format();
                    let (params, sketch): (SketchParams, LegacySketch) =
                        bincode::deserialize(bytes)?;
                    Ok((params, sketch.0))
                }
                SketchFormat::Unsupported(version) => Err(format!(
                    "Sketch database record {} was written by a newer version of skani \
                     (sketch format {}). Please upgrade skani or re-run `skani sketch`.",
                    index,
                    version.map_or("unknown".to_string(), |v| v.to_string())
                )
                .into()),
            }
        } else {
            Err(format!("Sketch index out of bounds: {}", index).into())
        }
    }

    /// Get the number of sketches in the database
    pub fn sketch_count(&self) -> usize {
        self.index.len()
    }

    /// Get the number of sketches in the database
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Check if the database is empty
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }
}

/// Check if a directory contains a consolidated sketch database
pub fn is_consolidated_db(database_dir: &str) -> bool {
    let pair_exists = |sketches: &str, index: &str| {
        Path::new(&format!("{}/{}", database_dir, sketches)).exists()
            && Path::new(&format!("{}/{}", database_dir, index)).exists()
    };
    pair_exists(SKETCHES_DB, INDEX_DB) || pair_exists(SKETCHES_DB_LEGACY, INDEX_DB_LEGACY)
}

/// Check if a directory contains separate sketch files (legacy format)
pub fn has_separate_sketches(database_dir: &str) -> bool {
    if let Ok(entries) = std::fs::read_dir(database_dir) {
        for entry in entries.flatten() {
            if let Some(file_name) = entry.file_name().to_str() {
                if file_name.ends_with(".sketch") {
                    return true;
                }
            }
        }
    }
    false
}
