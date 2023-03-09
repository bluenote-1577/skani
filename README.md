# skani - accurate, fast nucleotide identity calculation for MAGs and databases

## Introduction

**skani** is a program for calculating average nucleotide identity (ANI) from DNA sequences (contigs/MAGs/genomes) for ANI > ~80%.

skani uses an approximate mapping method without base-level alignment to get orthology to get ANI. It is magnitudes faster than BLAST based methods and almost as accurate. skani offers:

1. **Accurate ANI calculations for MAGs**. skani is accurate for incomplete and medium-quality metagenome-assembled genomes (MAGs). Sketching methods (e.g. Mash), which may underestimate ANI for incomplete MAGs.

2. **Aligned fraction results**. skani outputs the fraction of genome aligned, whereas sketching methods do not. 

3. **Fast computations**. Indexing/sketching is ~ 3x faster than Mash, and querying is about 25x faster than FastANI (but slower than Mash). 

4. **Efficient database search**. Querying a genome against a preprocessed database of >65000 prokaryotic genomes takes a few seconds with a single processor and ~6 GB of RAM. Constructing a database from genome sequences takes a few minutes to an hour.

##  Install

#### Option 1: Build from source

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) programming language and associated tools such as cargo are required and assumed to be in PATH.
2. A c compiler (e.g. GCC)
3. make

Building takes a few minutes (depending on # of cores).

```sh
git clone https://github.com/bluenote-1577/skani
cd skani

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo
skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta

# If ~/.cargo doesn't exist use below commands instead
#cargo build --release
#./target/release/skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta
```

#### Option 2: Pre-built x86-64 linux statically compiled executable

We offer a pre-built statically compiled executable for x86-64 linux systems. That is, if you're on a x86-64 linux system, you can just download the binary and run it without installing anything. 

For using the latest version of skani: 

```sh
wget https://github.com/bluenote-1577/skani/releases/download/latest/skani
chmod +x skani
./skani -h
```

Note: the binary is compiled with a different set of libraries (musl instead of glibc), possibly impacting performance (slightly). Probably not a huge deal. 

See the [Releases](https://github.com/bluenote-1577/skani/releases) page for obtaining specific versions of skani.


#### Option 3: Conda (conda version: 0.1.0 - source version: 0.1.0)

```sh
conda install -c bioconda skani
```

Note: skani is being developed quickly and the conda version may be outdated; see the version status above. 

## Quick start

```sh
# compare two genomes for ANI. 
# all options take -t for multi-threading.
skani dist genome1.fa genome2.fa -t 5

# compare multiple genomes
skani dist -q query1.fa query2.fa -r reference1.fa reference2.fa -o all-to-all_results.txt

# construct database and do memory-efficient search
skani sketch genomes_to_search/* -o database
skani search query1.fa query2.fa ... -d database

# use sketch from "skani sketch" output as drop-in replacement
skani dist database/query.fa.sketch database/ref.fa.sketch

# construct distance matrix for all genomes in folder
skani triangle genome_folder/* > skani_ani_matrix.txt

# we provide a script in this repository for clustering/visualizing distance matrices.
# requires python3, seaborn, scipy/numpy, and matplotlib.
python scripts/clustermap_triangle.py skani_ani_matrix.txt 

```

## Tutorials and manuals

### [skani basic usage information](https://github.com/bluenote-1577/skani/wiki/skani-basic-usage-guide)

For more information about using the specific skani subcommands, see the [guide linked above](https://github.com/bluenote-1577/skani/wiki/skani-basic-usage-guide). 

### skani tutorials

1. #### [Tutorial: setting up a 65000 prokaryotic genome database to search against](https://github.com/bluenote-1577/skani/wiki/Tutorial:-setting-up-a-65000-genome-database-to-search-against)
2. #### [Tutorial: strain-level clustering of MAGs using skani, and why Mash/FastANI have issues](https://github.com/bluenote-1577/skani/wiki/Tutorial:-strain-and-species-level-clustering-of-MAGs-with-skani-triangle)

### [skani advanced usage information](https://github.com/bluenote-1577/skani/wiki/skani-advanced-usage-guide)

See the advanced usage guide linked above for more information about topics such as:

* optimizing sensitivity/speed of skani
* using skani for long-reads
* making skani for memory efficient for huge data sets

## Output

If the resulting aligned fraction for the two genomes is < 15%, no output is given. 

**In practice, this means that only results with > ~82% ANI are reliably output** (with default parameters). See the [skani advanced usage guide](https://github.com/bluenote-1577/skani/wiki/skani-advanced-usage-guide) for information on how to compare lower ANI genomes. 

The default output for `search` and `dist` looks like
```
Ref_file	Query_file	ANI	Align_fraction_ref	Align_fraction_query	Ref_name	Query_name
refs/e.coli-EC590.fasta	refs/e.coli-K12.fasta	99.39	93.95	93.37	NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome	NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence
```
- Ref_file: the filename of the reference.
- Query_file: the filename of the query.
- ANI: the ANI.
- Aligned_fraction_query/reference: fraction of query/reference covered by alignments.
- Ref/Query_name: the id of the first record in the reference/query file.

## Citation

Jim Shaw and Yun William Yu. Fast and robust metagenomic sequence comparison through sparse chaining with skani. bioRxiv (2023).  https://doi.org/10.1101/2023.01.18.524587. Submitted.

##  Updates

### v0.1.0 released - 2023-02-07. 

We added new experiments on the revised version of our preprint (Extended Data Figs 11-14). We show skani has quite good AF correlation with MUMmer, and that it works decently on simple eukaryotic MAGs, especially with the `--slow` option (see below). 

#### Major

* **ANI debiasing added** - skani now uses a debiasing step with a regression model trained on MAGs to give more accurate ANIs. Old version gave robust, but slightly overestimated ANIs, especially around 95-97% range. Debiasing is enabled by default, but can be turned off with ``--no-learned-ani``.
* **More accurate aligned fraction** - chaining algorithm changed to give a more accurate aligned fraction (AF) estimate. The previous version had more variance and underestimated AF for certain assemblies.

#### Minor

* **Small contig/genome defaults made better** - should be more sensitive so that they don't get filtered by default.
* **Repetitive k-mer masking made better** - smarter settings and should work better for eukaryotic genomes; shouldn't affect prokaryotic genomes much.
* **`--fast` and `--slow` mode added** - alias for `-c 200` and `-c 30` respectively.
* **More non x86_64 builds should work** - there was a bug before where skani would be dysfunctional on non x86_64 architectures. It seems to at least build on ARM64 architectures successfully now.

## Feature requests, issues

skani is actively being developed by me ([Jim Shaw](https://jim-shaw-bluenote.github.io/)). I'm more than happy to accommodate simple feature requests (different types of outputs, etc). Feel free to open an issue with your feature request on the github repository. If you catch any bugs, please open an issue or e-mail me (e-mail on my website). 

## Calling skani from rust or python

### Rust API

If you're interested in using skani as a rust library, check out the minimal example here: https://github.com/bluenote-1577/skani-lib-example. The documentation is currently minimal (https://docs.rs/skani/0.1.0/skani/) and I guarantee no API stability. 

### Python bindings 

If you're interested in calling skani from python, see the [pyskani](https://github.com/althonos/pyskani) python interface and bindings to skani written by [Martin Larralde](https://github.com/althonos). Note: I am not personally involved in the pyskani project and do not offer guarantees on correctness of the outputs. 
