v0.1.3

- Fixed a bug where memory was blowing up in dist and triangle when the marker-index was activated.
- For all modes, implemented writing outputs during processing instead of storing all results until the end of the command. 
- Changed the marker index hash table population method. Used to overestimate memory usage slightly.
- New help message for marker parameters. Turns out that for small genomes, having more markers may make filtering significantly better. 
- Added -i option to sketch so you can sketch individual records in multifastas -- does not work for search yet though, only for sketching. 

v0.1.2

- Added medium preset.
- Added distance argument in triangle for distance instead of similarity matrices.
- Changed --marker-index option to --no-marker-index, which is a much more sane option. 


