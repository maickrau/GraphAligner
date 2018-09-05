This experiment downloads the E.Coli reference genome and builds four types of graphs from the first 10000 bp. Reads are simulated and aligned to all graphs with the bitvector and cell-by-cell algorithms, and the runtimes are reported.

1. The binaries for Snakemake (for running the experiment pipeline), BCALM v2.2 (for building the tangle graph) and PBSIM (for simulating the reads) must first be installed
   - https://snakemake.readthedocs.io/en/stable/
   - https://github.com/GATB/bcalm/releases
   - https://github.com/GATB/bcalm/blob/master/scripts/convertToGFA.py
   - https://github.com/pfaucon/PBSIM-PacBio-Simulator
2. Run `git clone git@github.com:maickrau/GraphAligner.git && cd GraphAligner && git checkout WabiExperiments` 
3. Run `make`
4. Edit `WabiExperimentSnake/config.yaml`: add the paths to bcalm and pbsim binaries, bcalm graph conversion script and the pbsim simulation parameter file
5. Run `cd WabiExperimentSnake && snakemake all`
6. The results will be in `WabiExperimentSnake/results/ref10000_summary.txt`
7. For multiple runs: copy `ref10000_summary.txt` somewhere, run `rm WabiExperimentSnake/results/ref10000_*.txt` and repeat steps 5-7

The files `ref10000_(graphname).txt` will have:
- Graph type (linear/tree/forest or DAG or cyclic)
- The size of the graph (nodes, edges): Nodes/edges is one-character nodes (the size reported in the paper), Collapsed nodes/edges is nodes with linear parts collaped into one node
- Total time to align all reads with bitvector and cellbycell, and the ratio
- Sanity check that the bitvector and cellbycell produce the same results

The output of `ref10000_summary.txt` will have:
- Total time to align all reads with bitvector and the cellbycell for each graph and the ratio for that run
- There is also one extra graph, the onechar graph, which is the linear graph with each node split into one character labels. This was not included in the paper because cellbycell has 64x node overhead compared to bitvector, so it's not a fair comparison. The linear graph is the fair version of this graph. The bitvector gets about a 16x speedup for the onechar graph compared to 10x for the linear graph.
