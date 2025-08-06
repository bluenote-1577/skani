use crate::params::*;
use crate::types::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
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
pub struct SketchDbWriter {
    concat_file: BufWriter<File>,
    index: Vec<IndexEntry>,
    current_offset: u64,
}

/// Reader for accessing consolidated sketch databases
pub struct SketchDbReader {
    mmap: memmap2::Mmap,
    index: HashMap<String, (u64, u64)>, // file_name -> (offset, length)
}

impl SketchDbWriter {
    /// Create a new consolidated sketch database writer
    pub fn new(output_dir: &str) -> Result<Self, Box<dyn std::error::Error>> {
                
        let concat_path = format!("{}/sketches.db", output_dir);
        let concat_file = BufWriter::new(File::create(concat_path)?);
        
        Ok(SketchDbWriter {
            concat_file,
            index: Vec::new(),
            current_offset: 0,
        })
    }

    /// Add a sketch to the consolidated database
    pub fn add_sketch(&mut self, sketch_params: &SketchParams, sketch: &Sketch) -> Result<(), Box<dyn std::error::Error>> {
        // Serialize the sketch data
        let serialized = bincode::serialize(&(sketch_params, sketch))?;
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
        let index_path = format!("{}/index.db", output_dir);
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
        let index_path = format!("{}/index.db", database_dir);
        let index_file = File::open(index_path)?;
        let index_reader = BufReader::new(index_file);
        let index_vec: Vec<IndexEntry> = bincode::deserialize_from(index_reader)?;

        // Convert to HashMap for O(1) lookups
        let mut index = HashMap::new();
        for entry in index_vec {
            index.insert(entry.file_name, (entry.offset, entry.length));
        }

        // Memory map the main database file
        let concat_path = format!("{}/sketches.db", database_dir);
        let concat_file = File::open(concat_path)?;
        let mmap = unsafe { memmap2::Mmap::map(&concat_file)? };

        info!("Loaded consolidated sketch database with {} sketches", index.len());
        Ok(SketchDbReader { mmap, index })
    }

    /// Get a sketch by file name
    pub fn get_sketch(&self, file_name: &str) -> Result<(SketchParams, Sketch), Box<dyn std::error::Error>> {
        use std::time::Instant;
        
        let total_start = Instant::now();
        
        if let Some(&(offset, length)) = self.index.get(file_name) {
            let start = offset as usize;
            let end = start + length as usize;
            
            // Time the mmap slice access (potential disk I/O)
            let bytes = &self.mmap[start..end];
            
            // Time the deserialization (memory copy + parsing)
            let deserialize_start = Instant::now();
            let (params, sketch) = bincode::deserialize(bytes)?;
            let deserialize_time = deserialize_start.elapsed();
            
            let total_time = total_start.elapsed();
                        
            Ok((params, sketch))
        } else {
            Err(format!("Sketch not found: {}", file_name).into())
        }
    }

    /// Get all sketch file names in the database
    pub fn get_sketch_names(&self) -> Vec<String> {
        self.index.keys().cloned().collect()
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
    let concat_path = format!("{}/sketches.db", database_dir);
    let index_path = format!("{}/index.db", database_dir);
    Path::new(&concat_path).exists() && Path::new(&index_path).exists()
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
