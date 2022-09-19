use crate::constants::*;
use crate::seeding;
use crate::types::*;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq::Reader as FastqReader;

pub fn fasta_to_sketch(ref_file: &str, k: usize, c: usize, use_syncs: bool) -> Sketch {
    let mut new_sketch = Sketch::default();
    new_sketch.k = k;
    new_sketch.file_name = ref_file.split('/').last().unwrap().to_string();
    let is_fastq = new_sketch.file_name.chars().last().unwrap() == 'q';
    dbg!(is_fastq);
    let mut i = 0;
    if is_fastq {
        let mut reader = FastqReader::from_file(ref_file).unwrap();
        let records = reader.records();
        for record in records{
//        while let Some(record) = records{
            let record = record.expect("Error reading record");
//            let contig = record.id().unwrap();
            let contig = record.id();
            new_sketch.contigs.push(contig.to_string());
            let seq = record.seq();
            let seq_len = seq.len();
            if seq_len < TRIANGLE_LENGTH_CUTOFF {
//                                println!("Sequence {} skipped because the length {} is less than {}", contig, seq.len(), TRIANGLE_LENGTH_CUTOFF);
                continue;
            }
            new_sketch.total_sequence_length += seq.len();
            if use_syncs {
                let adj_c =
                    seeding::open_sync_seeds(seq, k, c, i as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            } else {
                let adj_c = seeding::fmh_seeds(seq, k, c, i as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            }
            i += 1;
        }
    } else {
        let mut reader = FastaReader::from_file(ref_file).unwrap();
        let records = reader.records();
//        while let Some(record) = reader.next() {
        for record in records{
            let record = record.expect("Error reading record");
            let contig = record.id();
            new_sketch.contigs.push(contig.to_string());
            let seq = record.seq();
            let seq_len = seq.len();
            if seq_len < TRIANGLE_LENGTH_CUTOFF {
//                                println!("Sequence {} skipped because the length {} is less than {}", contig, seq.len(), TRIANGLE_LENGTH_CUTOFF);
                continue;
            }
            new_sketch.total_sequence_length += seq.len();
            if use_syncs {
                let adj_c =
                    seeding::open_sync_seeds(seq, k, c, i as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            } else {
                let adj_c = seeding::fmh_seeds(seq, k, c, i as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            }
            i += 1;
        }
    }
    return new_sketch;
}

pub fn fastq_to_multiple_sketch(
    ref_file: &str,
    k: usize,
    c: usize,
    use_syncs: bool,
) -> Vec<Sketch> {
    let mut new_sketches = vec![];
    let file_name = ref_file.split('/').last().unwrap().to_string();
    let is_fastq = file_name.chars().last().unwrap() == 'q';
    if is_fastq {
        let mut reader = FastqReader::from_file(ref_file).unwrap();
//        while let Some(record) = reader.next() {
        for record in reader.records(){
            let mut new_sketch = Sketch::default();
            new_sketch.k = k;
            let record = record.expect("Error reading record");
            let contig = record.id();
            new_sketch.contigs.push(contig.to_string());
            let seq = record.seq();
            new_sketch.total_sequence_length = seq.len();
            if use_syncs {
                let adj_c =
                    seeding::open_sync_seeds(seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            } else {
                let adj_c = seeding::fmh_seeds(seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            }
            new_sketches.push(new_sketch);
        }
    } else {
        let mut reader = FastaReader::from_file(ref_file).unwrap();
//        while let Some(record) = reader.next() {
        for record in reader.records(){
            let mut new_sketch = Sketch::default();
            new_sketch.k = k;
            let record = record.expect("Error reading record");
            let contig = record.id();
            new_sketch.contigs.push(contig.to_string());
            let seq = record.seq();
            new_sketch.total_sequence_length = seq.len();
            if use_syncs {
                let adj_c =
                    seeding::open_sync_seeds(seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            } else {
                let adj_c = seeding::fmh_seeds(seq, k, c, 0 as u32, &mut new_sketch.kmer_seeds);
                new_sketch.c_adj = adj_c;
            }

            new_sketches.push(new_sketch);
        }
    }
    return new_sketches;
}
