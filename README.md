# skani - accurate, robust, and fast nucleotide/amino acid identity calculation for MAGs and databases

## Introduction

**skani** is a software package for calculating average nucleotide identity (ANI) or average amino acid identity (AAI) for metagenomic data. skani is designed for pairs of genomes or MAGs (metagenome-assembled contigs) with > 85% ANI and > 60% AAI and total sequence length > 50kb. 

The main advantages of skani compared to other methods such as FastANI or Mash are

1. **Robustness to fragmentation and incompleteness** in MAGs for comparing similar species. This is done 
by the usage of sparse k-mer chaining to find approximate alignments. 

2. **Extremely fast**. Indexing/sketching is ~ 2.5x faster than Mash, and querying is about 20x faster than FastANI (but slower than Mash).

3. **Memory efficient**. skani can search a huge database of genomes by storing partial information and loading the full index only if a 
potential match is detected.


### Requirements

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.

### Install

```
git clone https://github.com/bluenote-1577/skani
cd skani
cargo build --release
./target/release/skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta
```

`cargo build --release` builds the **skani** binary, which is found in the ./target/release/ directory. 

## Using skani

### skani dist - simple ANI/AAI calculation

```
#simple ANI calculation
skani dist genome1.fa genome2.fa 

#AAI computation
skani dist -a genome1.fa genome2.fa 

#all-to-all, 20 threads
skani dist -q query1.fa query2.fa -r ref1.fa ref2.fa -t 20

#query each record in a multi-fasta (--qi for query, --ri for reference)
skani dist --qi -q query1.fa -r ref1.fa
```

`dist` computes ANI/AAI between all queries and all references. If the resulting aligned fraction for the two genomes (see Outputs) is < 15% for ANI or 5% for AAI, no output is given. With skani default parameters, this means that only computations
with > ~85% ANI and > ~60% AAI are output. 

`dist` loads all reference and query genomes into memory. 
If you're searching against a database, `search` can use much less memory (see below). 

### skani search - memory-efficient ANI/AAI database queries

```
# use -a while sketching for AAI computation instead
skani sketch genome1.fa genome2.fa ... -o sketch_folder 
skani search -d sketch_folder query1.fa query2.fa ...
```
`search` is a memory efficient method of calculating ANI/AAI against a large reference database. Searching against
the GTDB database (> 65000 genomes) takes only 4.5 GB of memory using `search`. This is achieved by only
fully loading genomes that pass a filter into memory, and discarding the index after each query is done. 

If for some reason you're querying many small sequences (> 1000 small sequences), the loading step will dominate the ANI comparison, so consider 
using `dist` instead. 

`search` requires all reference genomes to be sketched first using `skani sketch` and output into a new folder. **The parameters
for `search` are obtained from the parameters used for the `sketch` option**, so if you sketch for AAI using the `-a` option, you
can only use `search` for AAI. 

If querying many genomes (> 100) against a large database (> 10000 genomes) consider using the --fast-screen option. This takes longer on start-up, but speeds up subsequent queries. 

### skani triangle - all-to-all ANI/AAI compution 
```
skani triangle genome1.fa genome2.fa genome3.fa -o lower_triangle_matrix.txt 
skani triangle -l list_of_genomes.txt -o sparse_output --sparse 
```

`triangle` outputs a lower-triangular matrix in phyllip format. 
It avoids doing n^2 computations and only does n(n-1)/2 computations as opposed to `dist`. 
Use `--sparse` to output in the same format as `search` or `dist`. 

### skani sketch - storing sketches/indices on disk
```
skani sketch genome1.fa genome2.fa ... -o sketch_folder 
skani dist sketch_folder/genome1.fa.sketch sketch_folder/genome2.fa.sketch
```

`sketch` computes the sketch of a genome and stores it in a new folder. For each file `genome.fa`, two new files `genome.fa.sketch` and `genome.fa.markers`
are created in the output folder. The .sketch files can be used as drop-in substitutes for fasta files, and the .marker files are only used
in `search`. 

## Output formats

The default output for `search` and `dist` looks like
```
Ref_file	Query_file	ANI	Align_fraction_query	Align_fraction_reference	ANI_95_percentile	ANI_5_percentile	Ref_name	Query_name
data/e.coli-K12.fasta	data/e.coli-EC590.fasta	0.9939	0.9400	0.9342	0.9960	0.9919	NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence	NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome
```
- Ref_file: the filename of the reference.
- Query_file: the filename of the query.
- ANI/AAI: the ANI or AAI estimate.
- Aligned_fraction_query/reference: percentage of query/reference covered by alignments.
- ANI_95/5_percentile: heuristic 95% and 5% confidence intervals. In my experience, they are relatively accurate for ANI calculations between 95-99.9% but not for AAI. 
- Ref/Query_name: the id of the first contig in the reference/query file.

## Advanced

### Adjusting memory/speed tradeoffs 

### Different types of ANI/AAI estimation



