use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn test_sketch() {
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir1")
        .arg("-c")
        .arg("50")
        .arg("-k")
        .arg("13")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .arg("./test_files/e.coli-W.fasta.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir2")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir3")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir_aai")
        .arg("-a")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir1/")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .assert();
    assert.success().code(0);
}

#[test]
fn test_search() {
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);


    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir")
        .arg("-a")
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
        .arg("5").output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    println!("ANI search test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let ani = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[11]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[12]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 0.98);
    //assert!(af_q > 0.80);
    assert!(af_r > 0.80);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/o157_reads.fastq")
        .arg("--qi")
        .arg("--ql")
        .arg("./test_files/query_list.txt")
        .assert();
    assert.failure().code(2);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("--qi")
        .arg("--ql")
        .arg("./test_files/query_list.txt")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("search")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir/")
        .arg("./test_files/MN-03.fa")
        .arg("--marker-index")
        .arg("-n")
        .arg("5").output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    println!("AAI search test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    let ani = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[11]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[12]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 0.70);
    //assert!(af_q > 0.80);
    assert!(af_r > 0.50);
}
#[test]
fn test_dist() {
    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

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
    let aai = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[11]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[12]
        .parse::<f64>()
        .unwrap();
    assert!(aai > 0.98);
    assert!(af_q > 0.80);
    assert!(af_r > 0.80);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-EC590.fasta")
        .arg("./test_files/e.coli-K12.fasta")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[11]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[12]
        .parse::<f64>()
        .unwrap();
    assert!(ani > 0.992);
    assert!(ani < 1.00);
    assert!(af_q > 0.90);
    assert!(af_r > 0.90);

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
    let aai = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[11]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[12]
        .parse::<f64>()
        .unwrap();
    println!("AAI E.coli-klebsiella test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    assert!(aai > 0.70);
    assert!(aai < 0.85);
    assert!(af_q > 0.40);
    assert!(af_r > 0.40);
    assert!(af_q < 0.90);
    assert!(af_r < 0.90);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let out = cmd
        .arg("dist")
        .arg("./test_files/e.coli-K12.fasta")
        .arg("./test_files/MN-03.fa")
        .arg("-c")
        .arg("30")
        .output();
    let out_line = std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap();
    let ani = out_line.split('\t').collect::<Vec<&str>>()[10]
        .parse::<f64>()
        .unwrap();
    let af_q = out_line.split('\t').collect::<Vec<&str>>()[11]
        .parse::<f64>()
        .unwrap();
    let af_r = out_line.split('\t').collect::<Vec<&str>>()[12]
        .parse::<f64>()
        .unwrap();
    println!("ANI E.coli-klebsiella test");
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    assert!(aai > 0.70);
    assert!(aai < 0.85);
    assert!(af_q > 0.20);
    assert!(af_r > 0.20);
    assert!(af_q < 0.70);
    assert!(af_r < 0.70);

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
        .arg("-a")
        .assert();
    assert.success().code(0);

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

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-l")
        .arg("./test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir_aai/")
        .arg("-a")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let cmd = cmd
        .arg("dist")
        .arg("-r")
        .arg("./tests/results/test_sketch_dir1/e.coli-EC590.fasta.sketch")
        .arg("-q")
        .arg("./test_files/o157_reads.fastq")
        .arg("--qi")
        .arg("--robust")
        .arg("--marker-index");
    let std_bytes = &cmd.output().as_ref().unwrap().stdout.clone()[0..500];
    let ste_bytes = cmd.output().as_ref().unwrap().stderr.clone();
    let stdout = std::str::from_utf8(&std_bytes).unwrap();
    let stderr = std::str::from_utf8(&std_bytes).unwrap();
    println!("read test");
    println!("{}", stdout);
    assert!(stdout.split('\t').collect::<Vec<&str>>().len() > 10);
    cmd.assert().success().code(0);

    //println!("{}", std::str::from_utf8(&cmd.output().as_ref().unwrap().stderr).unwrap());
}

#[test]
fn test_triangle() {
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
    let assert = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
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
        .output();

    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stdout).unwrap()
    );
    println!(
        "{}",
        std::str::from_utf8(&out.as_ref().unwrap().stderr).unwrap()
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

    let mut cmd = Command::cargo_bin("skani").unwrap();
    let assert = cmd
        .arg("triangle")
        .arg("-l")
        .arg("./test_files/query_list.txt")
        .arg("-o")
        .arg("./tests/results/output")
        .arg("-E")
        .arg("--min-aligned-fraction")
        .arg("0.03")
        .assert();
    assert.success().code(0);
}
