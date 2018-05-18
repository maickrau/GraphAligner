1. Snakemake (for running the experiment pipeline), BCALM (for building the tangle graph) and PBSIM (for simulating the reads) must first be installed
https://snakemake.readthedocs.io/en/stable/
https://github.com/GATB/bcalm
https://github.com/pfaucon/PBSIM-PacBio-Simulator
2. git clone git@github.com:maickrau/GraphAligner.git && cd GraphAligner && git checkout WabiExperiments
3. mkdir obj && mkdir bin && make bin/Aligner
4. Get a reference genome in .fasta format for the experiments and put it in the WabiExperimentSnake directory
5. Edit WabiExperimentSnake/config.yaml: add the name of the reference genome, add the paths to bcalm and pbsim
6. cd WabiExperimentSnake && snakemake all
7. The results will be in WabiExperimentSnake/results/(genomename)_summary.txt where genomename is the file name of the reference genome
8. For multiple runs: copy (genomename)_summary.txt somewhere, rm WabiExperimentSnake/results/(genomename)_*.txt and repeat steps 6-8

The files (genomename)_(graphname).txt will have:
-Graph type (linear/tree/forest or DAG or cyclic)
-The size of the graph (nodes, edges): Nodes/edges is one-character nodes (the size reported in the paper), Collapsed nodes/edges is nodes with linear parts collaped into one node
-Time of the alignment for bitvector and cellbycell, and the ratio
-Sanity check that the bitvector and cellbycell produce the same results
The output of (genomename)_summary.txt will have:
-Times for the bitvector and the cellbycell for each graph and the ratio for that run
-There is also one extra graph, the onechar graph, which is the linear graph with each node split into one character labels. This was not included in the paper because cellbycell has 64x node overhead compared to bitvector, so it's not a fair comparison. The linear graph is the fair version of this graph. The bitvector gets about a 16x speedup for the onechar graph compared to 10x for the linear graph.
