### v0.2.0 released - 2023-09-26

#### BREAKING

* --learned-ani feature was buggy before and now removed. 

#### Major

* Major bug found: debiasing for ANI was turned off if there were > 5000 queries present in skani search and skani dist. This bug is fixed now. 

#### Minor

* The rust API is changing in this version. Not published to Cargo yet (waiting on https://github.com/DDOtten/partitions/pull/3 to be published to crates...)

### v0.1.5 released - 2023-09-01

#### Major

Improved "N" character support: 

* changed query-reference selection method slightly via a slight hack, using marker seeds to estimate reference length instead. This makes it so NNN characters are not counted. 
* Now seeds with "N" characters present are no longer indexed. 

#### Minor
* --robust now uses the learned ANI debiasing procedure by default. 

### v0.1.4 released - 2023-06-14

#### Major
* skani triangle had a bug where if more than 5000 queries were present and --sparse or -E was not specified, the intermediate batch of 5000 queries would be written in sparse mode. 
* skani triangle -o was giving different upper triangle matrix instead of lower triangle (skani triangle > res gives lower triangle). Matrices are consistently lower triangle now.
* Changed to lto = true for release mode. I see anywhere from a 5-10% speedup for this.

#### Minor
* Changed some dependencies so no more dependencies on old crates that will deprecate. 

### v0.1.3 released - 2023-05-09 

#### Major
* Fixed a bug where memory was blowing up in `dist` and `triangle` when the marker-index was activated. For big datasets, there could be > 100 GBs of wasted memory. 
* skani now outputs intermediate results after processing each batch of 5000 queries. **This will mean that outputs may no longer be deterministically ordered if there are > 5000 genomes**, but you can sort the output file to get deterministic outputs, i.e ``skani triangle *.fa | sort -k 3 -n > sorted_skani_result.txt`` will guarantee deterministic output order. 

#### Minor 
* Changed the marker index hash table population method. Used to overestimate memory usage slightly.
* New help message for marker parameters. Turns out that for small genomes, having more markers may make filtering significantly better. 
* Added -i option to sketch so you can sketch individual records in multifastas -- does not work for search yet though, only for sketching. 

### v0.1.2 released - 2023-04-28.

Small fixes.

* Added `--medium` pre-set, which is just `-c 70`. Seems to work okay for comparing fragmented genomes. 
* **BREAKING**: Changed `--marker-index` to `--no-marker-index` as a more sane option. 
* Added `--distance` option to `skani triangle` to output distance matrix (i.e. 100 - ANI) instead of similarity matrix. 
* Misc. help message fixes

### v0.1.1 released - 2023-04-09. 

Small fixes.

* Made aligned fraction in `triangle mode` a full matrix by default. This is not a symmetric matrix since AF is not symmetric. 
* Misc. help message fixes 

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
