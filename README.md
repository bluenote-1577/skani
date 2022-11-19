# skani - accurate, fast nucleotide/amino acid identity calculation for MAGs and databases

## Introduction

**skani** is a software package for calculating average nucleotide identity (ANI) or average amino acid identity (AAI) for metagenomic data. skani is designed for pairs of genomes or MAGs (metagenome-assembled contigs) with > 85% ANI and > 60% AAI and total sequence length > 20kb. 

skani uses an approximate alignment method without base-level alignment. It is magnitudes faster than BLAST based methods and almost as accurate. skani offers:

1. **Accurate ANI calculations for similar MAGs**. Other methods, such as Mash, give estimates are not accurate when MAGs are < 90% complete. 

2. **Extremely fast**. Indexing/sketching is ~ 2.5x faster than Mash, and querying is about 20x faster than FastANI (but slower than Mash). 

3. **Efficient database search**. Querying a genome against a pre-sketched GTDB database (>65000 genomes) for ANI takes a few seconds with a single processor and ~4.5 GB of RAM, almost as fast as Mash.

4. **Efficient AAI calculation**. AAI calculation is about 10-20x slower than ANI, but much faster than any other BLAST+Prodigal based method. The sketches for AAI can get quite large but can be stored on disk and queried efficiently.

### Requirements and Install

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.

```
git clone https://github.com/bluenote-1577/skani
cd skani
cargo build --release
./target/release/skani dist refs/e.coli-EC590.fasta refs/e.coli-K12.fasta
```

`cargo build --release` builds the **skani** binary, which is found in the ./target/release/ directory. 

## Using skani

### skani sketch - storing sketches/indices on disk
```
skani sketch genome1.fa genome2.fa ... -o sketch_folder 
skani sketch -l list_of_genomes.txt -o sketch_folder
skani sketch -a/--aai genome1.fa genome2.fa ... -o aai_sketch_folder

skani dist sketch_folder/genome1.fa.sketch sketch_folder/genome2.fa.sketch
```

`sketch` computes the sketch (a.k.a index) of a genome and stores it in a new folder. For each file `genome.fa`, two new files `genome.fa.sketch` and `genome.fa.markers`
are created in the output folder. 

The .sketch files can be used as drop-in substitutes for fasta files. 


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

`dist` computes ANI/AAI between all queries and all references. `dist` loads all reference and query genomes into memory. 
If you're searching against a database, `search` can use much less memory (see below). With default settings, the entire GTDB database (65000 bacterial genomes) takes about 95 GB of memory to store in memory. 

### skani search - memory-efficient ANI/AAI database queries

```
# use -a while sketching for AAI computation instead
skani sketch genome1.fa genome2.fa ... -o sketch_folder -t (threads)
skani search -d sketch_folder query1.fa query2.fa ...  -t (threads) -o output.txt
```
`search` is a memory efficient method of calculating ANI/AAI against a large reference database. Searching against
the GTDB database (> 65000 genomes) takes only 4.5 GB of memory using `search`. This is achieved by only
fully loading genomes that pass a filter into memory, and discarding the index after each query is done. 

If for some reason you're querying many small sequences (> 1000 small sequences), the loading step will dominate the ANI comparison, so consider 
using `dist` instead if you have enough RAM. 

`search` requires all reference genomes to be sketched first using `skani sketch` and output into a new folder. **The parameters
for `search` are obtained from the parameters used for the `sketch` option**, so if you sketch for AAI using the `-a` option, you
can only use `search` for AAI. 

### skani triangle - all-to-all ANI/AAI compution 
```
skani triangle genome1.fa genome2.fa genome3.fa -o lower_triangle_matrix.txt
skani triangle -l list_of_genomes.txt -o sparse_output --sparse 
```

`triangle` outputs a lower-triangular matrix in phyllip format. 
It avoids doing n^2 computations and only does n(n-1)/2 computations as opposed to `dist`. 
Use `--sparse` to output in the same format as `search` or `dist`. 

## Output

If the resulting aligned fraction for the two genomes is < 15% for ANI or 5% for AAI, no output is given. This can be changed, see the -h options.

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

For `triangle`, if the genomes fail to meet the alignment fraction cutoff, 0 is output. The default output for `triangle` is a lower-triangular ANI/AAI matrix to stdout, and a lower-triangular aligned-fraction matrix to the skani_matrix.af file. The name can be changed using the -o option for all outputs. 

## Advanced

### Adjusting memory/speed tradeoffs 

#### Marker index 
The ``--marker-index`` option is available for `skani dist` and `skani search`. This loads all marker k-mers into a hash table for constant time filtering. Building the table takes a bit of time, and the table itself is ~10 GB for 60000 genomes with default parameters. This is turned off if less than 100 query files are input, and turned on automatically otherwise. If you're querying less than 100 genomes (or if you're using the --qi option), you need to turn this on manually. 

#### Adjusting c
If you want skani to run faster, the main parameter to adjust is the `-c` parameter. skani's speed and memory efficiency is inversely proportional to c as it controls the sampling rate of the k-mer seeds. As a default, c = 120 for ANI and c = 15 for AAI. For genomes of ANI > 95%, c can be comfortably made higher up to 200. AAI TODO

However, decreasing `c` **may not necessarily improve ANI/AAI accuracy for > 85% ANI genomes** since default parameters are built around c = 120 and c = 15.

### Comparing low-ANI/AAI genomes. 

skani focuses on ANI/AAI comparisons for genomes with > 85% ANI and > 60% AAI. To get more accurate results for low ANI/AAI values, one should use a lower value for `c`. 

For example, the supplied genome `refs/MN-03.fa` is a Klebsiella Pneumoniae genome, and running ``skani dist refs/MN-03.fa refs/e.coli-K12.fa`` returns nothing, because the two genomes do not have a good enough alignment. However, ``skani dist refs/MN-03.fa refs/e.coli-K12.fa -c 30`` returns an ANI estimate of ~79%. 

The aligned fraction output is also less accurate on low-similarity genomes, so decreasing `c` will make it more accurate. 

### ANI calculations for small genomes/reads

skani is not necessarily designed for comparing long-reads or small contigs, but it seems to work relatively well for ANI when the reads/contigs are long enough. 

- skani can not classify short-reads. Use a taxonomic classifier such as kraken for this.
- skani can not compare sets of short-reads. Use Mash or sourmash for this.
- skani has not been tested for AAI on long-reads

For small contigs or long-reads, here are some suggestions:

1. **Make sure to use `--marker-index` option for comparing against large databases.**
2. Use the --qi option for `skani search` or `skani dist`.
3. `skani dist` will be much faster than `skani search`, since the bottleneck will be loading genomes into memory. 

For parameters:

1. The default marker size -m is set to 1000, so we take one marker per every 1000 k-mers. A good rule of thumb is that you want at least 30 markers on average, so set -m = avg_read_length / 30. 
2. Set -c to 60 for noisy-long reads (mean identity < 95), or keeping c = 120 for higher quality reads (or if memory is an issue). Longer + High identity => higher -c. 
3. skani currently loads the entire file into memory instead of processing one read at a time. Consider splitting large sets of reads. 

## Citation

TODO
