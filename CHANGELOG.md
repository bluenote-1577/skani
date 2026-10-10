### v0.4.0 - 2026-10-10 (Breaking changes; new database / sketches and 2x faster)

#### Major

* **Roughly 2x faster.** For all commands,`triangle / dist / search`. Results are unchanged (see the one exception under Minor).
* **~30% less memory.** 295 MB -> 210 MB on the benchmark above; 1105 MB -> 812 MB with `--slow`. A resident sketch index is about half its former size.
* **BREAKING: the sketch format changed.** Sketches written by 0.4 **cannot be read by skani 0.3.x or earlier.** The consolidated database files are named `sketches.db2` and `index.db2` so that an older binary fails to find them rather than misreading them; if you point skani 0.3.x at a 0.4 database it will exit with `Problem reading .../markers.bin` and produce no output.
* **Old sketches still work.** 0.4 reads pre-0.4 `.sketch` files, databases and `markers.bin`, and prints one warning per run. They are converted on every load, so re-running `skani sketch` is worth it: loading is substantially faster afterwards.

#### Minor

* Chain scores are now computed in integer arithmetic rather than floating point. This is exact, not an approximation, but it removes an unspecified tie-break: when two anchors tie for the best chain score, the winner used to depend on the internal iteration order of the union-find structure skani used. In testing this changed 1 output row in 4950 at `--medium`, by 0.01 ANI. Other modes were unaffected.
* Dropped the `partitions` dependency.

### v0.3.2 - 2026-05-12

### Major

- Commands with `-i` specified (triangle, sketch, etc) would not work on non-x86 platforms due to a logic error. This is now fixed.

### Minor

- UTF-8 indexing issue should be fixed when file paths have multibyte characters

### v0.3.1 - 2025-10-11

### Minor
* Fixed `--both-min-af` bug in `skani search`

### v0.3.0 released - 2025-08 (Breaking changes)

#### Major

* Engineering changes: skani takes ~30-40% less memory than before but is ~5-10% slower. Results should still be identical to before. 
* BREAKING: Changed sketching output. Now all sketches are concatenated into a searchable database file by default. Original behavior can be resotored via `--separate-files`. 
* New command line options `--both-min-af` and `--short-header` 

#### Minor

* Refactored the commandline backend. Small deviations in commandline behavior may be present, hopefully not big bugs were introduced. 

### v0.2.2 released - 2024-07-04

#### Major

* added the `--small-genomes` preset. This is just an alias for `-c 30 -m 200 --faster-small`. This makes skani much faster when comparing hundreds of thousands of small genomes. 

#### Minor

* fixed a bug where `skani triangle --full-matrix` gave different results between STDOUT and `-o` (thanks to Florian Plaza Onate)
* added a `--diagonal` option (suggested by Antonio Camargo) to print diagonal entries for sparse and lower-triangular distance matrices
* added a warning to use `--faster-small` when comparing too many contigs (e.g. viruses, plasmids). 

### v0.2.1 released - 2023-10-11

More consistent support for small contigs and sequences. 

#### Major

* --faster-small option included in dist and triangle. 

Genomes (and contigs with the --i, --ri, --qi options) with less than 20 marker k-mers are not screened according to the -s option. This makes skani for senative for small sequences, but can hamper performance on very large datasets with lots of small genomes/contigs. 

This heuristic can now be disabled with the `--faster-small` option. 

#### Minor

* skani's version is now displayed properly
* Added some error messages for degenerate cases (and more testing)
* We found that the statically built binary can be a lot slower in certain cases. File i/o may be an issue for the binary version. A note is now added in the README.

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
