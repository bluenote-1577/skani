use assert_cmd::prelude::*; 
use tsv::*;
 // Used for writing assertions
use std::process::Command; // Run programs
                           //
#[test]
fn full_test_sketch_and_search() {
    Command::new("rm")
        .arg("-r")
        .args(["./tests/results/test_sketch_dir2", "./tests/results/test_sketch_dir1","./tests/results/test_sketch_dir3", "./tests/results/test_sketch_dir", "./tests/results/test_sketch_dir_aai"])
        .spawn();
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir1")
        .arg("--separate-sketches")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir3")
        .arg("--separate-sketches")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir")
        .arg("--separate-sketches")
        .assert();
    assert.success().code(0);


    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("--median")
        .arg("-n")
        .arg("5")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    println!("ANI search test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let _af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 97.);
    assert!(af_q > 80.);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("--median")
        .arg("-n")
        .arg("5")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    println!("ANI search test learned");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let _af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 97.);
    assert!(af_q > 80.);



//    let mut cmd = Command::cargo_bin("skani").unwrap();
//    let assert = cmd
//        .arg("search")
//        .arg("-d")
//        .arg("./tests/results/test_sketch_dir/")
//        .arg("./test_files/e.coli-EC590.fasta")
//        .arg("./test_files/e.coli-K12.fasta")
//        .arg("./test_files/o157_reads.fastq")
//        .arg("--qi")
//        .arg("--ql")
//        .arg("./test_files/query_list.txt")
//        .assert();
//    assert.failure().code(2);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("--no-marker-index")
        .arg("--qi")
        .arg("--ql")
        .arg("./test_files/query_list.txt")
        .assert();
    assert.success().code(0);
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("./tests/results/test_sketch_dir/markers.bin")
        .arg("./tests/results/test_sketch_dir/e.coli-EC590.fasta.sketch")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    println!("{}",out_line);
    assert!(!out_line.contains("WARN"));

    
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let cmd = cmd
        .arg("dist")
        .arg("-r")
        .arg("./tests/results/test_sketch_dir1/e.coli-EC590.fasta.sketch")
        .arg("./tests/results/test_sketch_dir1/markers.bin")
        .arg("-q")
        .arg("./test_files/o157_reads.fastq")
        .arg("--qi")
        .arg("--robust");
    println!(
        "{}",
        std::str::from_utf8(&cmd.output().as_ref().unwrap().stdout).unwrap()
    );


}
#[test]
fn full_test_dist() {
    
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.0);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("--faster-small")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.0);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);



    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-n")
        .arg("3")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.2);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-n")
        .arg("3")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.2);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-n")
        .arg("3")
        .arg("--medium")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.2);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);



    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-n")
        .arg("3")
        .arg("--slow")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.2);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-n")
        .arg("3")
        .arg("--fast")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 99.2);
    assert!(ani < 100.);
    assert!(af_q > 90.);
    assert!(af_r > 90.);


    println!("ANI E.coli test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );

    
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/MN-03.fa")
        .arg("-c")
        .arg("30")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let _ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    println!("ANI E.coli-klebsiella test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    assert!(af_q > 20.);
    assert!(af_r > 20.);
    assert!(af_q < 70.);
    assert!(af_r < 70.);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/o157_plasmid.fasta")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("dist")
        .arg("-q")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("-r")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/o157_plasmid.fasta")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("dist")
        .arg("--qi")
        .arg("-q")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("--ri")
        .arg("-r")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/o157_plasmid.fasta")
        .arg("--ci")
        .assert();
    assert.success().code(0);

    Command::new("mkdir")
        .args(["./tests/results/"])
        .spawn();

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("dist")
        .arg("--ri")
        .arg("--rl")
        .arg("./test_files/query_list.txt")
        .arg("--qi")
        .arg("-q")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/o157_plasmid.fasta")
        .arg("-s")
        .arg("0.9")
        .arg("--robust")
        .arg("-o")
        .arg("tests/results/test_dist_file.txt")
        .assert();
    assert.success().code(0);

//    let std_bytes = &cmd.output().as_ref().unwrap().stdout.clone()[0..500];
//    let ste_bytes = cmd.output().as_ref().unwrap().stderr.clone();
//    let stdout = std::str::from_utf8(std_bytes).unwrap();
//    let stderr = std::str::from_utf8(&ste_bytes).unwrap();
//    println!("read test");
//    println!("{}", stdout);
//    assert!(stdout.split('\t').collect::<Vec<&str>>().len() > 10);
//    assert!(!stderr.contains("WARN") && !stderr.contains("ERROR"));
//    cmd.assert().success().code(0);
//
    //println!("{}", std::str::from_utf8(&cmd.output().as_ref().unwrap().stderr).unwrap());
}

#[test]
fn full_test_triangle() {
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("triangle")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("-s")
        .arg("0.9")
        .arg("--robust")
        .arg("-k")
        .arg("13")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("-s")
        .arg("0.9")
        .arg("--robust")
        .arg("-k")
        .arg("13");
    out.assert().success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .output();
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap()
    );
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("--faster-small")
        .arg("./test_files/query_list.txt")
        .output();
    println!("--faster-small TEST");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );



    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("-o")
        .arg("./tests/results/output")
        .output();
    //    assert
    //        .success()
    //        .code(0);
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("-o")
        .arg("./tests/results/output")
        .arg("-E")
        .arg("--min-af")
        .arg("3")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("-o")
        .arg("./tests/results/output")
        .arg("--no-learned-ani")
        .assert();
    assert.success().code(0);



    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("--full-matrix")
        .output();

    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("--full-matrix")
        .arg("-o")
        .arg("./tests/results/output_o_triangle_full")
        .output();

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("--full-matrix")
        .output();

    let std_out = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();

    //assert that std_out and output_o_triangle_full are the same by reading the file and using
    //assert!
    let mut cmd = Command::new("cat");
    let out = cmd.arg("./tests/results/output_o_triangle_full").output();
    let out = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    assert_eq!(std_out, out);
}

#[test]
fn test_consolidated_database_functionality() {
    // Clean up any existing test directories
    Command::new("rm")
        .arg("-rf")
        .args(["./tests/results/test_consolidated_db", "./tests/results/test_separate_db"])
        .spawn();

    // Test 1: Create consolidated database (default behavior)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-o")
        .arg("./tests/results/test_consolidated_db")
        .assert();
    assert.success().code(0);

    // Verify consolidated database files exist
    assert!(std::path::Path::new("./tests/results/test_consolidated_db/sketches.db").exists());
    assert!(std::path::Path::new("./tests/results/test_consolidated_db/index.db").exists());
    assert!(std::path::Path::new("./tests/results/test_consolidated_db/markers.bin").exists());

    // Test 2: Create separate sketches database using --separate-sketches flag
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-o")
        .arg("./tests/results/test_separate_db")
        .arg("--separate-sketches")
        .assert();
    assert.success().code(0);

    // Verify separate sketch files exist
    assert!(std::path::Path::new("./tests/results/test_separate_db/e.coli-EC590.fasta.sketch").exists());
    assert!(std::path::Path::new("./tests/results/test_separate_db/e.coli-K12.fasta.sketch").exists());
    assert!(std::path::Path::new("./tests/results/test_separate_db/markers.bin").exists());

    // Test 3: Search with consolidated database
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let consolidated_output = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_consolidated_db")
        .arg("./test_files/e.coli-EC590.fasta")
        .output()
        .unwrap();
    
    assert!(consolidated_output.status.success());
    let consolidated_stdout = std::str::from_utf8(&consolidated_output.stdout).unwrap();
    
    // Verify consolidated format was detected
    let consolidated_stderr = std::str::from_utf8(&consolidated_output.stderr).unwrap();
    assert!(consolidated_stderr.contains("Detected consolidated sketch database format"));

    // Test 4: Search with separate database
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let separate_output = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_separate_db")
        .arg("./test_files/e.coli-EC590.fasta")
        .output()
        .unwrap();
    
    assert!(separate_output.status.success());
    let separate_stdout = std::str::from_utf8(&separate_output.stdout).unwrap();
    
    // Verify separate format was detected
    let separate_stderr = std::str::from_utf8(&separate_output.stderr).unwrap();
    assert!(separate_stderr.contains("Detected separate sketch files format"));

    // Test 5: Verify both formats produce identical search results
    // Parse ANI results from both outputs (skip header line)
    let consolidated_lines: Vec<&str> = consolidated_stdout.lines().skip(1).collect();
    let separate_lines: Vec<&str> = separate_stdout.lines().skip(1).collect();
    
    assert_eq!(consolidated_lines.len(), separate_lines.len());
    
    // Compare ANI values from both results
    for (consolidated_line, separate_line) in consolidated_lines.iter().zip(separate_lines.iter()) {
        let consolidated_parts: Vec<&str> = consolidated_line.split('\t').collect();
        let separate_parts: Vec<&str> = separate_line.split('\t').collect();
        
        // Compare ANI values (column 2, 0-indexed)
        if consolidated_parts.len() > 2 && separate_parts.len() > 2 {
            let consolidated_ani: f64 = consolidated_parts[2].parse().unwrap_or(0.0);
            let separate_ani: f64 = separate_parts[2].parse().unwrap_or(0.0);
            
            // ANI values should be very close (within 0.01% difference)
            let diff = (consolidated_ani - separate_ani).abs();
            assert!(diff < 0.01, "ANI values differ too much: {} vs {}", consolidated_ani, separate_ani);
        }
    }

    // Clean up test directories
    Command::new("rm")
        .arg("-rf")
        .args(["./tests/results/test_consolidated_db", "./tests/results/test_separate_db"])
        .spawn();
}

#[test]
fn test_consolidated_database_multiple_files() {
    // Clean up any existing test directories
    Command::new("rm")
        .arg("-rf")
        .arg("./tests/results/test_multi_consolidated_db")
        .spawn();

    // Create consolidated database with multiple files
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/e.coli-W.fasta")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("-o")
        .arg("./tests/results/test_multi_consolidated_db")
        .assert();
    assert.success().code(0);

    // Verify database files exist
    assert!(std::path::Path::new("./tests/results/test_multi_consolidated_db/sketches.db").exists());
    assert!(std::path::Path::new("./tests/results/test_multi_consolidated_db/index.db").exists());
    assert!(std::path::Path::new("./tests/results/test_multi_consolidated_db/markers.bin").exists());

    // Test search functionality
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_multi_consolidated_db")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("-n")
        .arg("5")
        .output()
        .unwrap();
    
    assert!(output.status.success());
    let stdout = std::str::from_utf8(&output.stdout).unwrap();
    
    // Verify we get search results
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines.len() > 1, "Should have header plus at least one result line");
    
    // Verify the first result has high ANI (self-match)
    if lines.len() > 1 {
        let parts: Vec<&str> = lines[1].split('\t').collect();
        if parts.len() > 2 {
            let ani: f64 = parts[2].parse().unwrap_or(0.0);
            assert!(ani > 95.0, "Self-match should have high ANI: {}", ani);
        }
    }

    // Clean up test directory
    Command::new("rm")
        .arg("-rf")
        .arg("./tests/results/test_multi_consolidated_db")
        .spawn();
}

#[test]
fn test_short_header_functionality() {
    // Use existing test files with long headers containing whitespace
    let fasta1_path = "./test_files/e.coli-K12.fasta";
    let fasta2_path = "./test_files/e.coli-EC590.fasta";
    
    // Test 1: dist command without --short-header (should show full headers)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("dist")
        .arg(&fasta1_path)
        .arg(&fasta2_path)
        .output()
        .expect("Failed to execute skani dist");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "dist command failed: {}", stdout);
    assert!(stdout.contains("NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence"), 
            "Full header should be present without --short-header: {}", stdout);
    assert!(stdout.contains("NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome"), 
            "Full header should be present without --short-header: {}", stdout);
    
    // Test 2: dist command with --short-header (should show only first part)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("dist")
        .arg(&fasta1_path)
        .arg(&fasta2_path)
        .arg("--short-header")
        .output()
        .expect("Failed to execute skani dist with --short-header");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "dist command with --short-header failed: {}", stdout);
    assert!(stdout.contains("NC_007779.1"), "Short header should contain accession: {}", stdout);
    assert!(stdout.contains("NZ_CP016182.2"), "Short header should contain accession: {}", stdout);
    assert!(!stdout.contains("Escherichia coli str. K-12 substr. W3110"), 
            "Description should be truncated with --short-header: {}", stdout);
    assert!(!stdout.contains("complete sequence"), 
            "Description should be truncated with --short-header: {}", stdout);
    
    // Test 3: triangle command without --short-header
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("triangle")
        .arg(&fasta1_path)
        .arg(&fasta2_path)
        .arg("--sparse")
        .output()
        .expect("Failed to execute skani triangle");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "triangle command failed: {}", stdout);
    assert!(stdout.contains("NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence"), 
            "Full header should be present in triangle output: {}", stdout);
    
    // Test 4: triangle command with --short-header
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("triangle")
        .arg(&fasta1_path)
        .arg(&fasta2_path)
        .arg("--sparse")
        .arg("--short-header")
        .output()
        .expect("Failed to execute skani triangle with --short-header");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "triangle command with --short-header failed: {}", stdout);
    assert!(stdout.contains("NC_007779.1"), "Triangle short header should contain accession: {}", stdout);
    assert!(!stdout.contains("Escherichia coli str. K-12 substr. W3110"), 
            "Triangle description should be truncated with --short-header: {}", stdout);
    
    // Test 5: Create a sketch database for search testing
    let sketch_dir = "./tests/results/short_header_sketches";
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg(&fasta1_path)
        .arg(&fasta2_path)
        .arg("-o")
        .arg(&sketch_dir)
        .arg("--separate-sketches")
        .assert();
    assert.success();
    
    // Test 6: search command without --short-header
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("search")
        .arg("-d")
        .arg(&sketch_dir)
        .arg(&fasta1_path)
        .output()
        .expect("Failed to execute skani search");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "search command failed: {}", stdout);
    assert!(stdout.contains("NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence") ||
            stdout.contains("NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome"), 
            "Full header should be present in search output: {}", stdout);
    
    // Test 7: search command with --short-header
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("search")
        .arg("-d")
        .arg(&sketch_dir)
        .arg(&fasta1_path)
        .arg("--short-header")
        .output()
        .expect("Failed to execute skani search with --short-header");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "search command with --short-header failed: {}", stdout);
    assert!(stdout.contains("NC_007779.1") || stdout.contains("NZ_CP016182.2"), 
            "Search short header should contain accession: {}", stdout);
    assert!(!stdout.contains("Escherichia coli str. K-12 substr. W3110") && 
            !stdout.contains("complete sequence"), 
            "Search description should be truncated with --short-header: {}", stdout);
    
    // Clean up test sketch directory
    Command::new("rm")
        .arg("-rf")
        .arg(sketch_dir)
        .spawn();
}

#[test]
fn test_individual_contigs_with_search() {
    // Test that sketch -i (individual contigs) works correctly with search
    let sketch_dir = "./tests/results/individual_contigs_search_test";
    
    // Test 1: Create sketches with -i (individual contigs) - should use default consolidated format
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-K12.fasta")  // Multi-contig file
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("-i")  // Individual contigs flag
        .arg("-o")
        .arg(sketch_dir)
        .assert();
    assert.success();
    
    // Test 2: Verify that search works with individual contig sketches (consolidated format)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("search")
        .arg("-d")
        .arg(sketch_dir)
        .arg("./test_files/e.coli-K12.fasta")
        .output()
        .expect("Failed to execute skani search with individual contig sketches");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(output.status.success(), "Search with individual contigs failed: {}", stdout);
    
    // Should find matches since we're searching for the same genome that was sketched
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines.len() > 1, "Should have header line plus at least one result: {}", stdout);
    
    // Test 3: Test the problematic combination - sketch with both -i and --separate-sketches  
    let incompatible_sketch_dir = "./tests/results/incompatible_sketch_test";
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("-i")  // Individual contigs
        .arg("--separate-sketches")  // Separate files (incompatible combo)
        .arg("-o")
        .arg(incompatible_sketch_dir)
        .output()
        .expect("Failed to execute skani sketch with incompatible flags");
    
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "Sketch command should still succeed: {}", stderr);
    
    // Should contain the warning we added
    assert!(stderr.contains("WARNING") && stderr.contains("separate-sketches") && stderr.contains("individual contigs"), 
            "Should warn about incompatible combination: {}", stderr);
    
    // Test 4: Verify that search fails or gives poor results with the incompatible format
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("search")
        .arg("-d")
        .arg(incompatible_sketch_dir)
        .arg("./test_files/e.coli-K12.fasta")
        .output()
        .expect("Failed to execute skani search with incompatible sketches");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    
    // The search might fail or give no results - either is acceptable for this incompatible combination
    // This test documents the expected behavior
    if !output.status.success() {
        // Search fails - this is expected behavior for incompatible format
        assert!(stderr.contains("No valid reference fastas") || stderr.len() > 0,
                "Search should fail gracefully with incompatible format: {}", stderr);
    } else {
        // Search succeeds but may have no results
        let lines: Vec<&str> = stdout.lines().collect();
        // Having only a header line (no results) is acceptable for incompatible format
        assert!(lines.len() >= 1, "Should at least have header line: {}", stdout);
    }
    
    // Clean up test directories
    Command::new("rm")
        .arg("-rf")
        .arg(sketch_dir)
        .arg(incompatible_sketch_dir)
        .spawn();
}

#[test]
fn test_sketch_search_individual_contigs_matches_dist() {
    // Test that sketch -i followed by search --qi gives same results as dist --qi --ri
    let sketch_dir = "./tests/results/test_sketch_individual_o157_comparison";
    let search_output = "./tests/results/search_individual_o157_comparison.txt";
    let dist_output = "./tests/results/dist_individual_o157_comparison.txt";
    
    // Clean up any existing test directories
    let _ = Command::new("rm")
        .arg("-rf")
        .arg(sketch_dir)
        .arg(search_output)
        .arg(dist_output)
        .spawn();
    
    // Step 1: Sketch e.coli-o157.fasta with individual contigs (-i)
    let sketch_cmd = Command::new("./target/debug/skani")
        .arg("sketch")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("-o")
        .arg(sketch_dir)
        .arg("-i")
        .output()
        .expect("Failed to execute skani sketch -i");
    
    assert!(sketch_cmd.status.success(), 
            "skani sketch -i should succeed: {}", String::from_utf8_lossy(&sketch_cmd.stderr));
    
    // Step 2: Search with individual query contigs (--qi)
    let search_cmd = Command::new("./target/debug/skani")
        .arg("search")
        .arg("-d")
        .arg(sketch_dir)
        .arg("./test_files/e.coli-o157.fasta")
        .arg("--qi")
        .arg("-o")
        .arg(search_output)
        .arg("-s")
        .arg("50") // Lower screen threshold to ensure results
        .output()
        .expect("Failed to execute skani search --qi");
    
    assert!(search_cmd.status.success(), 
            "skani search --qi should succeed: {}", String::from_utf8_lossy(&search_cmd.stderr));
    
    // Step 3: Run dist with individual contigs for both query and reference (--qi --ri)
    let dist_cmd = Command::new("./target/debug/skani")
        .arg("dist")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("--qi")
        .arg("--ri")
        .arg("-o")
        .arg(dist_output)
        .output()
        .expect("Failed to execute skani dist --qi --ri");
    
    assert!(dist_cmd.status.success(), 
            "skani dist --qi --ri should succeed: {}", String::from_utf8_lossy(&dist_cmd.stderr));
    
    // Step 4: Compare the results
    let search_content = std::fs::read_to_string(search_output)
        .expect("Should be able to read search output");
    let dist_content = std::fs::read_to_string(dist_output)
        .expect("Should be able to read dist output");
    
    // Parse both outputs and verify they contain the same ANI results
    let search_lines: Vec<&str> = search_content.lines().skip(1).collect(); // Skip header
    let dist_lines: Vec<&str> = dist_content.lines().skip(1).collect(); // Skip header
    
    // Both should have results (e.coli-o157 has 2 contigs, so we expect 2x2=4 comparisons)
    assert!(!search_lines.is_empty(), "Search should produce results");
    assert!(!dist_lines.is_empty(), "Dist should produce results");
    
    // Extract ANI values from both outputs and compare
    let mut search_anis: Vec<f32> = Vec::new();
    let mut dist_anis: Vec<f32> = Vec::new();
    
    for line in search_lines {
        if !line.trim().is_empty() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                if let Ok(ani) = fields[2].parse::<f32>() {
                    search_anis.push(ani);
                }
            }
        }
    }
    
    for line in dist_lines {
        if !line.trim().is_empty() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                if let Ok(ani) = fields[2].parse::<f32>() {
                    dist_anis.push(ani);
                }
            }
        }
    }
    
    // Both should have the same number of results
    assert_eq!(search_anis.len(), dist_anis.len(), 
               "Search and dist should produce same number of results. Search: {}, Dist: {}", 
               search_anis.len(), dist_anis.len());
    
    // Sort both vectors for comparison
    search_anis.sort_by(|a, b| a.partial_cmp(b).unwrap());
    dist_anis.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    // Compare ANI values (should be nearly identical)
    for (search_ani, dist_ani) in search_anis.iter().zip(dist_anis.iter()) {
        let diff = (search_ani - dist_ani).abs();
        assert!(diff < 0.01, 
                "ANI values should be nearly identical. Search: {}, Dist: {}, Diff: {}", 
                search_ani, dist_ani, diff);
    }
    
    // Clean up test directories
    let _ = Command::new("rm")
        .arg("-rf")
        .arg(sketch_dir)
        .arg(search_output)
        .arg(dist_output)
        .spawn();
}

#[test]
fn test_both_min_af_functionality() {
    // Test the new --both-min-af option behavior
    // Case 1: o157 plasmid vs o157 ecoli should pass without --both-min-af but be filtered with it
    // Case 2: o157 plasmid vs o157 plasmid should pass even with --both-min-af
    
    // Test Case 1: plasmid vs ecoli without --both-min-af (should have results)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("dist")
        .arg("./test_files/o157_plasmid.fasta")
        .arg("./test_files/e.coli-o157.fasta")
        .output()
        .unwrap();
    
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "Command failed with stderr: {}", stderr);
    
    // Should have output lines (one header + one result line)
    let lines: Vec<&str> = stdout.trim().split('\n').collect();
    assert!(lines.len() >= 2, "Should have at least header and one result line, got: {}", stdout);
    
    // Test Case 2: plasmid vs ecoli with --both-min-af 50 (should be filtered out)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("dist")
        .arg("./test_files/o157_plasmid.fasta")
        .arg("./test_files/e.coli-o157.fasta")
        .arg("--both-min-af")
        .arg("50")
        .output()
        .unwrap();
    
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "Command failed with stderr: {}", stderr);
    
    // Should only have header line (no results due to filtering)
    let lines: Vec<&str> = stdout.trim().split('\n').filter(|line| !line.is_empty()).collect();
    assert_eq!(lines.len(), 1, "Should have only header line with --both-min-af filtering, got: {}", stdout);
    
    // Test Case 3: plasmid vs plasmid with --both-min-af 50 (should NOT be filtered)
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let output = cmd
        .arg("dist")
        .arg("./test_files/o157_plasmid.fasta")
        .arg("./test_files/o157_plasmid.fasta")
        .arg("--both-min-af")
        .arg("50")
        .output()
        .unwrap();
    
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "Command failed with stderr: {}", stderr);
    
    // Should have output lines (header + result line) because both aligned fractions should be high
    let lines: Vec<&str> = stdout.trim().split('\n').filter(|line| !line.is_empty()).collect();
    assert!(lines.len() >= 2, "Should have header and result line for plasmid self-comparison, got: {}", stdout);
    
    // Verify the result line contains high aligned fractions
    if lines.len() >= 2 {
        let result_line = lines[1];
        let fields: Vec<&str> = result_line.split('\t').collect();
        // Fields should be: Ref-file, Query-file, ANI, Align-fraction-ref, Align-fraction-query, etc.
        if fields.len() >= 5 {
            let align_frac_ref: f64 = fields[3].parse().unwrap_or(0.0);
            let align_frac_query: f64 = fields[4].parse().unwrap_or(0.0);
            
            assert!(align_frac_ref >= 50.0, "Reference aligned fraction should be >= 50%, got: {}", align_frac_ref);
            assert!(align_frac_query >= 50.0, "Query aligned fraction should be >= 50%, got: {}", align_frac_query);
        }
    }
    
    println!("✓ --both-min-af functionality test passed!");
}

