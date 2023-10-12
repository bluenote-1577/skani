use assert_cmd::prelude::*; 
use reflection::Reflection;
use std::collections::HashMap;
use tsv::*;
 // Used for writing assertions
use serial_test::serial;
use std::process::Command; // Run programs
                           //
#[derive(Default, Clone, Debug, PartialEq)]
pub struct AniResult{
    pub ani: f32,
    pub align_fraction_query: f32,
    pub align_fraction_ref: f32,
    pub ref_file: String,
    pub query_file: String,
    pub query_contig: String,
    pub ref_contig: String,
}

fn run_skani<'a>(args: &'a [&str], stderr: bool) -> String{
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .args(args)
        .output();
    let out_line;
    if stderr{
        out_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap().to_string();
    }
    else{
        out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap().to_string();
    }
    return out_line
}

fn get_result_from_out(tsv_res : &str) -> Vec<AniResult>{
    let lines : Vec<&str> = tsv_res.lines().collect();
    let mut ret = vec![];
    for line in &lines[1..]{
        let line = line.split('\t').collect::<Vec<&str>>();
        if line.len() < 7{
            return ret;
        }
        ret.push(AniResult{
            ani: line[2].parse::<f32>().unwrap(),
            align_fraction_query: line[4].parse::<f32>().unwrap(),
            align_fraction_ref: line[3].parse::<f32>().unwrap(),
            ref_file: line[0].to_string(),
            query_file: line[1].to_string(),
            query_contig: line[5].to_string(),
            ref_contig: line[6].to_string(),
        });
    }
    ret
}

#[test]
fn fast_test_small_genomes() {
    let out_line = run_skani(&["triangle", "-i", "-E", "./test_files/viruses.fna"], false);
    let results = get_result_from_out(&out_line);
    for result in results{
        assert!(result.ani > 99.0);
        assert!(result.ani < 99.9);
    }

    let out_line = run_skani(&["triangle", "-i", "-E", "./test_files/viruses.fna", "--faster-small"], false);
    let results = get_result_from_out(&out_line);
    for result in results{
        assert!(result.ani > 99.0);
        assert!(result.ani < 99.9);
    }


    let out_line = run_skani(&["triangle", "-i", "-E", "./test_files/o157_reads.fastq"], false);
    let results = get_result_from_out(&out_line);
    println!("o157_reads, number of results without faster-small: {}", results.len());

    let out_line = run_skani(&["triangle", "-i", "-E", "--faster-small", "./test_files/o157_reads.fastq"], false);
    let results = get_result_from_out(&out_line);
    println!("o157_reads, number of results WITH faster-small: {}", results.len());
}

#[test]
fn fast_test_screen(){

    //Screening doesn't matter when -m is too high for markers to exist. But we should still get
    //results
    let out_line = run_skani(&["triangle", "-i", "-E", "./test_files/o157_reads.fastq", "-s","95", "-m", "1000000"], false);
    let results_s95 = get_result_from_out(&out_line);
    println!("rescue small heuristic with screen {}", results_s95.len());

    let out_line = run_skani(&["triangle", "-i", "-E", "./test_files/o157_reads.fastq", "-m", "10000000"], false);
    let results_nos = get_result_from_out(&out_line);
    println!("rescue small heuristic no screen {}", results_nos.len());

    let out_line = run_skani(&["dist", "--qi",  "./test_files/o157_reads.fastq", "./test_files/e.coli-EC590.fasta", "--faster-small"], false);
    let results = get_result_from_out(&out_line);
    println!("dist no rescue, no screen {}", results.len());

    let out_line = run_skani(&["dist", "--qi",  "./test_files/o157_reads.fastq", "./test_files/e.coli-EC590.fasta"], false);
    let results = get_result_from_out(&out_line);
    println!("dist rescue, no screen {}", results.len());

    let out_line = run_skani(&["dist", "--qi",  "./test_files/o157_reads.fastq", "./test_files/e.coli-EC590.fasta", "-s", "95", "--faster-small"], false);
    let results_s95_faster = get_result_from_out(&out_line);
    println!("no rescue, screen {}", results_s95_faster.len());

    let out_line = run_skani(&["dist", "--qi",  "./test_files/o157_reads.fastq", "./test_files/e.coli-EC590.fasta", "--faster-small", "-m",  "10000000000"], false);
    let results_s95_faster_large = get_result_from_out(&out_line);
    println!("no rescue, screen, 0 markers {}", results_s95_faster_large.len());

    assert!(results_s95.len() == results_nos.len());
    assert!(results_s95.len() > 0);
    assert!(results_s95_faster_large.len() == 0);
}

#[test]
fn fast_show_degenerate_inputs(){
    let out_line = run_skani(&["dist", "file_that_doesn't_exist", "file_that_doesn't_exist"], true);
    assert!(out_line.contains("WARN"));
    println!("{}", out_line);

    let out_line = run_skani(&["dist", "./test_files/query_list.txt", "./test_files/viruses.fna"], true);
    assert!(out_line.contains("WARN"));
    println!("{}", out_line);

    let out_line = run_skani(&["dist", "./test_files/query_list.txt", "list.txt"], true);
    assert!(out_line.contains("WARN"));
    println!("{}", out_line);

    let out_line = run_skani(&["dist", "./test_files/empty_fasta.fa", "./test_files/empty_fasta.fa"], true);
    assert!(out_line.contains("WARN"));
    println!("{}", out_line);

    let out_line = run_skani(&["dist", "--qi", "--ri", "./test_files/empty_fasta.fa", "./test_files/empty_fasta.fa"], true);
    assert!(out_line.contains("WARN"));
    println!("{}", out_line);

    //Ns don't get sketched. 
    let out_line = run_skani(&["dist", "./test_files/all_ns.fa", "./test_files/all_ns.fa"], true);
    println!("{}", out_line);
    let results = get_result_from_out(&out_line);
    assert!(results.len() == 0);

}
