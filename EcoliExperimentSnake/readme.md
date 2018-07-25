This experiment downloads Illumina and Pacbio reads for E. Coli, builds a de Bruijn graph of the Illumina reads, selects Pacbio reads longer than 1000bp, randomly downsamples them into ~1.5x coverage and aligns them to the de Bruijn graph with the bitvector and cell-by-cell algorithms. Note: this experiment will take a long time to run since the entire DP matrix is calculated!

1. The binaries for Snakemake (for running the experiment pipeline), BCALM (for building the graph) and sra-toolkit (for downloading the Pacbio reads) must first be installed
   - https://snakemake.readthedocs.io/en/stable/
   - https://github.com/GATB/bcalm
   - https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
2. Run `git clone git@github.com:maickrau/GraphAligner.git && cd GraphAligner && git checkout WabiExperiments` 
3. Run `mkdir obj && mkdir bin && make bin/Aligner`
4. Edit `EcoliExperimentSnake/config.yaml`: add the paths to the bcalm binary and graph conversion script. You can also optionally change the k-mer size and minimum abundance for the de Bruijn graph.
5. Run `cd EcoliExperimentSnake && snakemake all`
6. The results will be in `EcoliExperimentSnake/results/results/k(k)_cov(cov)_summary.txt`

For multiple runs, it's a better idea to first run `snakemake tmp/pacbios_downsampled.fa tmp/graph_k(k)_cov(cov).gfa` and then start multiple parallel runs of `../bin/Aligner -g tmp/graph_k(k)_cov(cov).gfa -f tmp/pacbios_downsampled.fa > results.txt`. Replace (k) with the k-mer size of the graph and (cov) with the k-mer solidity threshold.

The output file will have the total time to align all reads for the bitvector and cellbycell, and the ratio.
