Binsanity2 is the upcoming new version of Binsanity designed to alleviate some of the memory constraints of the original iteration. 

The `Binsanity2-Beta` script has been tested and shown to produce similarly high quality results compared to original versions.

The average user should use the `--skip-kmeans` flag if you working with a small assembly that is < 25,000 contigs. 

If you have ample RAM (e.g at least 400-500GB) you can use the `--skip-kmeans` for anything <50,000 contigs.

This flag bypasses the Kmeans clustering step and tells the program to ONLY use Affinity Propagation. 

Above 50,000 contigs you will want the k-means clustering step to limit memory usage. 

There are two workflow bash scripts included here with the beta version. These bash scripts are updated version of the multi-pass binning procedure
initially implemented in Tully et al. 2018 (https://doi.org/10.1038/sdata.2017.203) for the Tara Oceans metagenomes.   

One is designed for large assemblies (>25,000 contigs) and one designed for small assemblies <25,000. The difference being that 
"binsanity2-beta-wf-smallassemblies.sh" workflow uses the `--skip-kmeans` flag and the `binsanity2-beta-wf-largeassemblies.sh` uses custom pre-set `-C` values. 
