use bio::io::fasta;
use crate::seeding;
use crate::types::*;

pub fn fasta_to_sketch(ref_file: &str, k: usize, c: usize) -> Sketch{
    let reader = fasta::Reader::from_file(ref_file).unwrap();
    let mut new_sketch = Sketch::default();
    new_sketch.file_name = ref_file.split('/').last().unwrap().to_string();
    for (i,result) in reader.records().enumerate(){
        let record = result.unwrap();
        let contig = record.id();
        new_sketch.contigs.push(contig.to_string());
        let seq = record.seq();
        new_sketch.total_sequence_length += seq.len();
        let adj_c = seeding::fmh_seeds(seq, k, c, i as u32, &mut new_sketch.kmer_seeds);
        new_sketch.c_adj = adj_c;
    }

    return new_sketch;
}
