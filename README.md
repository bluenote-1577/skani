# skani - accurate, fast nucleotide/amino acid identity calculation for MAGs and databases

## Introduction

**skani** is a program for calculating average nucleotide identity (ANI) or average amino acid identity (AAI) from microbial DNA sequences (contigs/MAGs/genomes). 

skani is designed species or genus-level ANI calculations, whereas the AAI mode offers more sensitivity. 

skani uses an approximate mapping method to get orthology without base-level alignment, and then estimates ANI/AAI. It is magnitudes faster than BLAST based methods and almost as accurate. skani offers:

1. **Accurate ANI calculations for MAGs**. skani is more accurate than approximate sketching methods, such as Mash. skani remains accurate for incomplete and medium-quality metagenome-assembled genomes (MAGs), whereas Mash loses accuracy when genomes are incomplete. 

2. **Fast computations**. Indexing/sketching is ~ 2.5x faster than Mash, and querying is about 25x faster than FastANI (but slower than Mash). 

3. **Efficient database search**. Querying a genome against a preprocessed GTDB database (>65000 genomes) takes a few seconds with a single processor and ~4.5 GB of RAM. Constructing a database from genome sequences takes only a few minutes. 

4. **Efficient AAI calculation on DNA sequences**. skani can calculate AAI between two DNA sequences (i.e. no gene prediction needed). AAI calculation for two genomes takes at most 1 second. Querying against a database can take a few minutes.

##  Install

#### Option 1: Build from source

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) programming language and associated tools such as cargo are required and assumed to be in PATH.

Building takes a few minutes (depending on # of cores).

```sh
git clone https://github.com/bluenote-1577/skani
cd skani
# make sure ~/.cargo exists
cargo install --path . --root ~/.cargo
skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta
```

#### If ~/.cargo is not present (non-standard rust installs)
```
cargo build --release
./target/release/skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta
```

#### Option 2: Pre-built linux binary

We offer a pre-built binary for 64-bit linux systems. This is convenient but may not be as up-to-date as the source. 

See the [Releases](https://github.com/bluenote-1577/skani/releases) page. 

```sh
wget https://github.com/bluenote-1577/skani/releases/download/binary-release/skani-linux-1.0.0-alpha.tar
tar -xf skani-linux-1.0.0-alpha.tar
cd skani-linux-1.0.0-alpha
./skani -h
```



## Quick start

All skani modes take the argument `-t` as number of threads (default: 3).

```sh
# compare two genomes for ANI. 
# all options take -t for multi-threading.
skani dist genome1.fa genome2.fa -t 5

# use -a or --aai for all commands to calculate AAI instead of ANI
skani dist -a genome1.fa genome2.fa > aai_results.txt

# compare multiple genomes
skani dist -q query1.fa query2.fa -r reference1.fa reference2.fa -o all-to-all_results.txt

# construct database and do memory-efficient search
skani sketch genomes_to_search/* -o database
skani search query1.fa query2.fa ... -d database

# construct distance matrix for all genomes in folder
skani triangle genome_folder/* > distance_matrix.txt

# we provide a script in this repository for clustering/visualizing distance matrices.
# requires python3, seaborn, scipy/numpy, and matplotlib.
python scripts/clustermap_triangle.py distance_matrix.txt 

```

### [skani commands usage information](https://github.com/bluenote-1577/skani/wiki/skani-basic-usage-guide)

See [skani commands usage information](https://github.com/bluenote-1577/skani/wiki/skani-basic-usage-guide) for an explanation of each of skani's subcommands, and which one to use for your situation.

### skani tutorials

1. #### [Tutorial: setting up a 65000 bacterial genome database to search against](https://github.com/bluenote-1577/skani/wiki/Tutorial:-setting-up-a-65000-genome-database-to-search-against)

### [skani advanced usage information](https://github.com/bluenote-1577/skani/wiki/skani-advanced-usage-guide)

For more information about topics such as:

* using skani for long-reads
* making skani for memory efficient for huge data sets
* optimizing skani for speed/memory tradeoffs

see [skani advanced usage information](https://github.com/bluenote-1577/skani/wiki/skani-advanced-usage-guide).

## Output

If the resulting aligned fraction for the two genomes is < 15% for ANI or 5% for AAI, no output is given. This can be changed, see the `--min-aligned-fraction` option.

**In practice, this means that only genomes with > ~82% ANI and > ~60% AAI are output** with default parameters. 

The default output for `search` and `dist` looks like
```
Ref_file	Query_file	ANI	Align_fraction_query	Align_fraction_reference	ANI_95_percentile	ANI_5_percentile	Ref_name	Query_name
data/e.coli-K12.fasta	data/e.coli-EC590.fasta	0.9939	0.9400	0.9342	0.9960	0.9919	NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence	NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome
```
- Ref_file: the filename of the reference.
- Query_file: the filename of the query.
- ANI/AAI: the ANI or AAI estimate.
- Aligned_fraction_query/reference: fraction of query/reference covered by alignments.
- ANI_95/5_percentile: heuristic 95% and 5% confidence intervals. IThey are relatively accurate for ANI calculations between 95-99.9% on prokaryotic MAGs/genomes, but not for AAI or small genomes. 
- Ref/Query_name: the id of the first contig in the reference/query file.

## Citation

Jim Shaw and Yun William Yu. Fast and robust metagenomic sequence comparison through sparse chaining with skani (2022). Submitted.
