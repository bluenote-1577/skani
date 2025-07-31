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
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir_aai")
        .arg("-a")
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
    let out = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir_aai/")
        .arg("./test_files/MN-03.fa")
        .arg("-n")
        .arg("5")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

    println!("AAI search test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let ani = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let _af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 70.);
    //assert!(af_q > 0.80);
    assert!(af_r > 30.);
    let err_line = std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap();
    assert!(!err_line.contains("WARN") && !err_line.contains("ERROR"));

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
        .arg("-a")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    println!("AAI E.coli test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let aai = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    assert!(aai > 98.);
    assert!(af_q > 80.);
    assert!(af_r > 80.);

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
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/MN-03.fa")
        .arg("-a")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    println!("{:?}", out_line.split('\t').collect::<Vec<&str>>());

    let aai = out_line.split('\t').collect::<Vec<&str>>()[8]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[9]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    println!("AAI E.coli-klebsiella test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    assert!(aai > 70.);
    assert!(aai < 85.);
    assert!(af_q > 30.);
    assert!(af_r > 30.);
    assert!(af_q < 90.);
    assert!(af_r < 90.);

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
    assert!(aai > 70.);
    assert!(aai < 85.);
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
    assert!(std::path::Path::new("./tests/results/test_consolidated_db/test_consolidated_db.db").exists());
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
    assert!(std::path::Path::new("./tests/results/test_multi_consolidated_db/test_multi_consolidated_db.db").exists());
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

