# GraphAligner

Seed-and-extend program for aligning long error-prone reads to genome graphs. For simple usage, see ["Running"](#running). For using the snakemake file, see ["Running the snakemake pipeline"](#running-the-snakemake-pipeline). For a description of the bitvector alignment extension algorithm, see https://www.biorxiv.org/content/early/2018/05/15/323063

### Installation

- Install miniconda https://conda.io/docs/user-guide/install/index.html
- `git clone https://github.com/maickrau/GraphAligner.git`
- `git submodule update --init --recursive`
- `conda env create -f CondaEnvironment.yml`
- `source activate GraphAligner`
- `make bin/Aligner`
- Optional: for the snakemake pipeline `make all`

### Running

Quickstart: `bin/Aligner -g graph_file -f read_file -a output_file.gam`

See [Parameters](#parameters), the option `bin/Aligner --help` and the subsections below for more information and options

#### File formats

The aligner's file formats are interoperable with the [vg toolkit](https://github.com/vgteam/vg/)'s file formats. Graphs can be inputed either in [.gfa format](https://github.com/GFA-spec/GFA-spec) or [.vg format](https://github.com/vgteam/vg/blob/master/src/vg.proto). Reads are inputed as .fasta or .fastq, either gzipped or uncompressed. Alignments are outputed in [.gam format](https://github.com/vgteam/vg/blob/master/src/vg.proto). Seeds can be inputed in [.gam format](https://github.com/vgteam/vg/blob/master/src/vg.proto).

#### Seed hits

The aligner has two built-in methods for finding seed hits: maximal unique matches (MUMs) (default) and maximal exact matches (MEMs). These modes use [MUMmer4](https://github.com/mummer4/mummer) to find matches between the read and nodes. Only matches entirely within a node are found. Use the parameter `--seeds-mum-count n` to use the `n` longest MUMs as seeds (or -1 for all MUMs), and `--seeds-mem-count n` for the `n` longest MEMs (or -1 for all MEMs). Use `--seeds-mxm-length n` to only use matches at least `n` characters long. If you are aligning multiple files to the same graph, use `--seeds-mxm-cache-prefix file_name_prefix` to store the MUM/MEM index to disk for reuse instead of rebuilding it each time.

Alternatively you can use any method to find seed hits and then import the seeds in [.gam format](https://github.com/vgteam/vg/blob/master/src/vg.proto) with the parameter `-s seedfile.gam`. The seeds must be passed as an alignment message, with `path.mapping[0].position` describing the position in the graph, `name` the name of the read and `query_position` the position in the read. Match length (`path.mapping[0].edit[0]`) is only used to order the seeds, with longer matches tried before shorter matches.

Alternatively you can use the parameter `--seeds-first-full-rows` to use the dynamic programming alignment algorithm on the entire first row instead of using seeded alignment. This is very slow except on tiny graphs, and not recommended.

#### Extension

The aligner uses a bitvector banded DP alignment algorithm to extend the seed hits. The DP matrix is calculated inside a certain area (the band), which depends on the extension parameters. Note that "bandwidth" in graph alignment does NOT directly correspond to bandwidth in linear alignment. The bandwidth parameter describes the maximum allowed score difference between the minimum score in a row and a cell, with cells whose score is higher than that falling outside the band. Generally the bandwidth parameters should be between 1-35.

The algorithm starts using the initial bandwidth. Should it detect that the alignment is incorrect, it will rewind and rerun with the ramp bandwidth parameter, aligning high-error parts of the read without slowing down alignment in low-error parts. The tangle effort parameter determines how much time the aligner spends inside complex cyclic subgraphs. If the size of the band grows beyond the tangle effort parameter, the aligner will use the current best alignment for the aligned prefix and move forward along the read. This might miss the optimal alignment.

### Running the snakemake pipeline

- Copy the folder Snakefiles/align_with_mummer somewhere, call it run_folder
- Acquire a graph and save it in run_folder/input/
- Acquire reads and save them in run_folder/input/
- Edit run_folder/config.yaml
  - Add the path of the repository folder
  - Add the input file names
  - Optionally edit the aligner parameters
- Run `source activate GraphAligner`
- Run `snakemake --cores [number of threads] all` in run_folder
- The output will be in run_folder/output/

#### Output

A summary of the alignment is in `output/stats_[readfile]_[graphfile].txt`. Individual alignments are in the files `output/aln_[readfile]_[graphfile]_all.gam`, `output/aln_[readfile]_[graphfile]_selected.gam` and `output/aln_[readfile]_[graphfile]_full_length.gam`. To view the `.gam` files, use [vg](https://github.com/vgteam/vg/). 

The file `_all.gam` contains all found alignments, usually including multiple incompatible alignments per read, eg. one base pair is aligned to multiple parts of the graph. The file `_selected.gam` is a set of non-overlapping alignments, so each read base pair is aligned to at most one location in the graph. The file `_full_length.gam` contains alignments where the entire read was aligned to one location.

### Parameters

- `-b` alignment bandwidth. Unlike linear alignment, this is the score difference between the minimum score in a row and the score where a cell falls out of the band. Values should be between 1-35.
- `-B` extra bandwidth. If a read cannot be aligned with the aligner bandwidth, switch to extra bandwidth at the problematic location. Values should be between 1-35.
- `-C` tangle effort. Determines how much effort the aligner spends on tangled areas. Higher values use more CPU and memory and have a higher chance of aligning through tangles. Lower values are faster but might return an inoptimal or a partial alignment. Use for complex graphs (eg. de Bruijn graphs of mammalian genomes) to limit the runtime in difficult areas. Values should be between 1'000 - 500'000.
- `--sloppy-optimizations` unsafe heuristics. Speeds up alignment in some cases but might lead to inoptimal results.
- `--all-alignments` output all alignments. Normally only a set of non-overlapping partial alignments is returned. Use this to also include partial alignments which overlap each others.
- `--high-memory` high memory mode. Runs a bit faster but uses a lot more memory

Suggested example parameters:
- Variation graph: `-b 35`
- RNA splice graph: `-b 35`
- Human de Bruijn graph: `-b 5 -B 10 -C 1000 --sloppy-optimizations`
- Bacterial de Bruijn graph: `-b 10 -B 20`

The parameters below are only relevant if manually running GraphAligner. If you are using the snakefile, you shouldn't do anything with these

- `-g` input graph. Format .gfa / .vg
- `-f` input reads. Format .fasta / .fastq / .fasta.gz / .fastq.gz
- `-s` input seeds. Format .gam
- `-t` number of threads
- `-a` output file name. Format .gam
- `--seeds-mum-count` MUM seeds. Use the longest n maximal unique matches for alignment. -1 for all MUMs
- `--seeds-mem-count` MEM seeds. Use the longest n maximal exact matches for alignment. -1 for all MEMs
- `--seeds-mxm-length` MUM/MEM minimum length. Don't use MUMs/MEMs shorter than n
- `--seeds-mxm-cache-prefix` MUM/MEM file cache prefix. Store the MUM/MEM index into disk for reuse. Recommended unless you are sure you won't align to the same graph multiple times
- `--seeds-first-full-rows` Don't use seeds. Instead use the DP alignment on the first row. The runtime depends on the size of the graph so this is very slow. Not recommended
