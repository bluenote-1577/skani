# skani - accurate, fast nucleotide/amino acid identity calculation for MAGs and databases

## Introduction

**skani** is a software package for calculating average nucleotide identity (ANI) or average amino acid identity (AAI) from microbial DNA sequences. skani is designed for species or genus-level ANI calculations, whereas the AAI mode offers more sensitivity.

skani uses an approximate mapping method to get orthology without base-level alignment, and then estimates ANI/AAI. It is magnitudes faster than BLAST based methods and almost as accurate. skani offers:

1. **Accurate ANI calculations for MAGs**. skani is more accurate than approximate sketching methods, such as Mash. skani remains accurate for incomplete and medium-quality metagenome-assembled genomes (MAGs), whereas Mash loses accuracy when genomes are incomplete. 

2. **Fast computations**. Indexing/sketching is ~ 2.5x faster than Mash, and querying is about 25x faster than FastANI (but slower than Mash). 

3. **Efficient database search**. Querying a genome against a preprocessed GTDB database (>65000 genomes) takes a few seconds with a single processor and ~4.5 GB of RAM, almost as fast as Mash. Constructing a database from genome sequences takes only a few minutes. 

4. **Efficient AAI calculation on DNA sequences**. skani can calculate AAI between two DNA sequences without requiring the user to predict genes. AAI calculation for two genomes takes at most 1 second. Querying against a database can take a few minutes.

### Requirements and Install

#### Option 1: Pre-built binary

We offer a pre-built binary for 64-bit linux systems. 

```sh
# TODO ENSURE THIS LINK WORKS ON RELEASE
wget https://github.com/bluenote-1577/skani/releases/download/binary-release/skani-linux-1.0.0-alpha.tar
tar -xf skani-linux-1.0.0-alpha.tar
cd skani-linux-1.0.0-alpha
./skani -h
```

#### Option 2: Build from source

1. [rust](https://www.rust-lang.org/tools/install) programming language and associated tools such as cargo are required and assumed to be in PATH.

```sh
git clone https://github.com/bluenote-1577/skani
cd skani
cargo build --release

# skani binary found in ./target/release
./target/release/skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta

# OPTIONAL: soft link to skani in PATH. 
#~/.cargo/bin is present for default rust installations. Link to another directory in PATH if needed 
ln -s  $PWD/target/release/skani ~/.cargo/bin/ 
skani -h
```

#### Quick start

```sh
# compare two genomes for ANI
skani dist genome1.fa genome2.fa

# compare two genomes for AAI
skani dist -a genome1.fa genome2.fa > aai_results.txt

# compare multiple genomes
skani dist -q query1.fa query2.fa -r reference1.fa reference2.fa -o all-to-all_results.txt

# construct database and do memory-efficient search
skani sketch genomes_to_search/* -o database
skani search query1.fa query2.fa ... -d database

# construct distance matrix 
skani triangle genome1.fa genome2.fa ... > distance_matrix.txt
# we provide a script in this repository for visualizing distance matrices.
# requires python3, seaborn, scipy/numpy, and matplotlib.
python scripts/clustermap_triangle.py distance_matrix.txt 
```

## Using skani

### skani sketch - storing sketches/indices on disk
```sh
# sketch genomes, output in sketch_folder
skani sketch genome1.fa genome2.fa ... -o sketch_folder 

# sketch a list of genomes specified in a file
skani sketch -l list_of_genomes.txt -o sketch_folder

# sketch for AAI instead of ANI (ANI by default).
skani sketch -a/--aai genome1.fa genome2.fa ... -o aai_sketch_folder

# use sketch file for computation
skani dist sketch_folder/genome1.fa.sketch sketch_folder/genome2.fa.sketch
```

`sketch` computes the sketch (a.k.a index) of a genome and stores it in a new folder. For each file `genome.fa`,  the new file `sketch_folder/genome.fa.sketch` is created. 

The `.sketch` files can be used as faster drop-in substitutes for fasta files. 

A special file `markers.bin` is also constructed and used specifically for the `search` command. 


### skani dist - simple ANI/AAI calculation

```sh
# all-to-all, 20 threads
skani dist -q query1.fa query2.fa -r ref1.fa ref2.fa -t 20

# query each record in a multi-fasta (--qi for query, --ri for reference)
skani dist --qi -q query1.fa -r ref1.fa
```

`dist` computes ANI/AAI between all queries and all references. `dist` loads all reference and query genomes into memory. 
If you're searching against a database, `search` can use much less memory (see below). With default settings, the entire GTDB database (65000 bacterial genomes) takes about 95 GB of memory to store in memory. 

### skani search - memory-efficient ANI/AAI database queries

```sh
# use -a while sketching for AAI computation instead
skani sketch genome1.fa genome2.fa ... -o sketch_folder -t (threads)

# query query1.fa, query2.fa, ... against sketches in sketch_folder
skani search -d sketch_folder query1.fa query2.fa ...  -t (threads) -o output.txt
```
`search` is a memory efficient method of calculating ANI/AAI against a large reference database. Searching against
the GTDB database (> 65000 genomes) takes only 4.5 GB of memory using `search`. This is achieved by only
fully loading genomes that pass a filter into memory, and discarding the index after each query is done. 

`search` requires all reference genomes to be sketched first using `skani sketch` and output into a new folder. **The parameters
for `search` are obtained from the parameters used for the `sketch` option**, so if you sketch for AAI using the `-a` option, you
can only use `search` for AAI. 

If you're querying many sequences, the file I/O step will dominate the running time, so consider 
using `dist` instead if you have enough RAM. 

### skani triangle - all-to-all ANI/AAI compution 
```sh
# all-to-all ANI comparison in lower-triangular matrix
skani triangle genome1.fa genome2.fa genome3.fa -o lower_triangle_matrix.txt

# output sparse matrix a.k.a an edge list of comparisons
skani triangle -l list_of_genomes.txt -o sparse_matrix.txt --sparse 

# output square matrix
skani triangle genome1.fa genom2.fa genome3.fa --full-matrix 
```

`triangle` outputs a lower-triangular matrix in [phyllip format](https://mothur.org/wiki/phylip-formatted_distance_matrix/). 
It avoids doing n^2 computations and only does n(n-1)/2 computations as opposed to `dist`. 

`triangle` loads all genome indices into memory. For doing comparisons on massive data sets, see the Advanced section for suggestions on reducing memory cost.

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

## Advanced

## Adjusting memory/speed tradeoffs 

### Marker index 
The ``--marker-index`` option is available for `skani dist` and `skani search`. This loads all marker k-mers into a hash table for constant time filtering. Building the table takes a bit of time, and the table itself is ~10 GB for 60000 genomes with default parameters. This is turned off if less than 100 query files are input or when using the `--qi` option. Otherwise, it is turned on automatically. 

### Adjusting c
If you want skani to run faster, the main parameter to adjust is the `-c` parameter. skani's speed and memory efficiency is inversely proportional to c, so increasing c by 2x means 2x faster and less memory. As a default, c = 120 for ANI and c = 15 for AAI. 

**ANI**: for genomes of ANI > 95%, c can be comfortably made higher up to 200 or even 300, but aligned fraction may get less accurate. 

**AAI**: for genomes of AAI > 65%, c can be made up to 30 and still relatively accurate. Results degrade a bit after c gets past 40. 

However, decreasing `c` **may not necessarily improve ANI/AAI accuracy for > 85% ANI genomes** since many other default algorithm parameters are designed these default values. Furthermore, increasing c means that distant genomes will no longer be comparable; see the section on **Comparing lower ANI/AAI genomes**. 

## All-to-all comparisons on massive data sets

`skani triangle` should be used for all-to-all comparisons on reasonably sized data sets. However, it loads all genome indices into memory, so RAM may be an issue. If RAM is an issue, consider: 
1. Pre-sketch using `skani sketch -l list_of_genomes.txt -o sketched_genomes` and run `skani search -d sketched_genomes -l list_of_genomes -o output` to do slower but low-memory all-to-all comparisons.
2. Consider using a lower marker density, `-m`. It defaults to 1000 but 2000 is reasonable for most bacterial genomes. 
3. Raising the `-c` parameter can help, see the above section on the `-c` parameter. 

## Comparing lower ANI/AAI genomes. 

skani focuses on ANI/AAI comparisons for genomes with > 85% ANI and > 60% AAI. To get more accurate results for low ANI/AAI values, one should use a lower value for `c`. 

For example, the supplied genome `refs/MN-03.fa` is a Klebsiella Pneumoniae genome, and running ``skani dist refs/MN-03.fa refs/e.coli-K12.fa`` returns nothing because the two genomes do not have a good enough alignment. However, ``skani dist refs/MN-03.fa refs/e.coli-K12.fa -c 30`` returns an ANI estimate of ~79%. 

For distant genomes, the aligned fraction output becomes more accurate as `c` gets smaller. However, decreasing `c` may *not necessarily* make high ANI calculations more accurate. Nevertheless, I would not recommend ANI comparisons for genomes with < 75% ANI, and advise using skani's AAI method instead, which is tuned for sensitive comparisons by default.


## ANI calculations for small genomes/reads

skani is not necessarily designed for comparing long-reads or small contigs, but it seems to work relatively well for ANI when the reads/contigs are long enough. 

- skani can not classify short-reads. Use a taxonomic classifier such as kraken for this.
- skani can not compare *collections* of short-reads. Use Mash or sourmash for this.
- skani can not compute AAI for long-reads.

For small contigs or long-reads, here are some suggestions:

1. Make sure to use the --qi option for `skani search` or `skani dist` if your contigs/reads are all in one file. 
2. `skani dist` will be much faster than `skani search`, since the bottleneck will be loading genomes into memory. 

For parameters:

1. The default marker size -m is set to 1000, so we take one marker per every 1000 k-mers. A good rule of thumb is that you want at least 20 markers on average, so set -m > avg_read_length / 20. 
2. Set -c to 60 for noisy-long reads (mean identity < 95), or keeping c = 120 for higher quality reads (or if memory is an issue). Longer + High identity => higher -c. 
3. skani currently loads the entire file into memory instead of processing one read at a time. Consider splitting large sets of reads. 

## Citation

Jim Shaw and Yun William Yu. Fast and robust metagenomic sequence comparison through sparse chaining with skani (2022). Submitted.
