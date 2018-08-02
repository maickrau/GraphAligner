# GraphAligner

Program for aligning long error-prone reads to genome graphs. For simple usage, see ["Running the alignment pipeline"](#running-the-alignment-pipeline). For a description of the bitvector algorithm, see https://www.biorxiv.org/content/early/2018/05/15/323063

### Installation

- Compile the aligner
  - Install zstr development libraries https://github.com/mateidavid/zstr
  - Install concurrentqueue development libraries https://github.com/cameron314/concurrentqueue
  - Install protobuf v3.0.0 development libraries https://github.com/google/protobuf/releases/tag/v3.0.0
  - Install sparsehash development libraries https://github.com/sparsehash/sparsehash
  - `make all`
- Install programs used by the snakemake pipeline
  - Install MUMmerv4 binaries https://github.com/mummer4/mummer
  - Install snakemake binaries https://snakemake.readthedocs.io/en/stable/

### Running the alignment pipeline

- Copy the folder Snakefiles/align_with_mummer somewhere, call it run_folder
- Acquire a graph and save it in run_folder/input/
- Acquire reads and save them in run_folder/input/
- Edit run_folder/config.yaml
  - Point the paths at the repository folder, vg and mummerv4
  - Add the input file names
  - Optionally edit the aligner parameters
- Run `snakemake --cores [number of threads] all` in run_folder
- The output will be in run_folder/output/

### Output

A summary of the alignment is in `output/stats_[readfile]_[graphfile].txt`. Individual alignments are in the files `output/aln_[readfile]_[graphfile]_all.gam`, `output/aln_[readfile]_[graphfile]_selected.gam` and `output/aln_[readfile]_[graphfile]_full_length.gam`. To view the `.gam` files, use [vg](https://github.com/vgteam/vg/). 

The file `_all.gam` contains all found alignments, usually including multiple incompatible alignments per read, eg. one base pair is aligned to multiple parts of the graph. The file `_selected.gam` is a set of non-overlapping alignments, so each read base pair is aligned to at most one location in the graph. The file `_full_length.gam` contains alignments where the entire read was aligned to one location.

### Parameters

- `-b` alignment bandwidth. Unlike linear alignment, this is the score difference between the minimum score in a row and the score where a cell falls out of the band. Values should be between 1-35.
- `-B` extra bandwidth. If a read cannot be aligned with the aligner bandwidth, switch to extra bandwidth at the problematic location. Values should be between 1-35.
- `-C` bandwidth limit. If the bandwidth grows beyond this value, stop aligning. Use for complex graphs (eg. de Bruijn graphs of mammalian genomes) to skip unalignable areas. Values should be between 100'000 - 10'000'000.
- `-l` low memory mode. Uses a lot less memory but runs a bit slower.

Suggested example parameters:
- Variation graph: `-b 35`
- RNA splice graph: `-b 35`
- Human de Bruijn graph: `-b 5 -B 10 -C 1000000`
- Bacterial de Bruijn graph: `-b 10 -B 20`

The parameters below are only relevant if manually running GraphAligner. If you are using the snakefile, you shouldn't do anything with these

- `-g` input graph. Format .gfa / .vg
- `-f` input reads. Format .fasta / .fastq / .fasta.gz / .fastq.gz
- `-s` input seeds. Format .gam
- `-t` number of threads
- `-a` output file name. Format .gam
