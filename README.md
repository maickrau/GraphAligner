# GraphAligner

Seed-and-extend program for aligning long error-prone reads to genome graphs. For a description of the bitvector alignment extension algorithm, see https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz162/5372677

### Installation

Install via [bioconda](https://bioconda.github.io/):

- Install miniconda https://conda.io/projects/conda/en/latest/user-guide/install/index.html
- `conda install -c bioconda graphaligner`

#### Compilation

Bioconda is the recommended installation method. If you however want to compile the aligner yourself, run these:

- Install miniconda https://conda.io/projects/conda/en/latest/user-guide/install/index.html
- `git clone https://github.com/maickrau/GraphAligner.git`
- `cd GraphAligner`
- `git submodule update --init --recursive`
- `conda env create -f CondaEnvironment.yml`
- `source activate GraphAligner`
- `make bin/GraphAligner`

Note that miniconda is only required during compilation and not during runtime. After compilation you can run the binary without the miniconda environment or copy the binary elsewhere.

If you want to compile without miniconda, you will need to install [boost](https://www.boost.org/), [mummer](https://github.com/mummer4/mummer), [protobuf and protoc](https://developers.google.com/protocol-buffers), [sdsl](https://github.com/simongog/sdsl-lite), [jemalloc](https://github.com/jemalloc/jemalloc) and [sparsehash](https://github.com/sparsehash/sparsehash).

### Running

Quickstart: `GraphAligner -g test/graph.gfa -f test/read.fa -a test/aln.gaf -x vg`

See [Parameters](#parameters), the option `GraphAligner --help` and the subsections below for more information and options

#### Test example

The command above outputs the following alignment in `test/aln.gaf`:

`read	71	0	71	+	>1>2>4	87	3	73	66	72	255	NM:i:6	dv:f:0.0833333	id:f:0.916667	cg:Z:4M1M2M1M1I37M1D5M1I5M1M13M`

which aligned the read to the nodes 1,2,4 with an identity of 91.66%. See [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) for more information about the output format. Alternatively try `-a aln.gam` for output compatible with [vg](https://github.com/vgteam/vg/).

The parameter `-x vg` uses a parameter preset for aligning reads to a variation graph. Other options are `-x dbg` for aligning to a de Bruijn graph.

#### File formats

The aligner's file formats are interoperable with [vg](https://github.com/vgteam/vg/)'s file formats. Graphs can be inputed either in [.gfa format](https://github.com/GFA-spec/GFA-spec) or [.vg format](https://github.com/vgteam/libvgio/blob/master/deps/vg.proto). Reads are inputed as .fasta or .fastq, either gzipped or uncompressed. Alignments are outputed in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) or [vg's alignment format](https://github.com/vgteam/libvgio/blob/master/deps/vg.proto), either as a binary .gam or JSON depending on the file name. Custom seeds can be inputed in [.gam format](https://github.com/vgteam/libvgio/blob/master/deps/vg.proto).

#### Seed hits

The aligner has three built-in methods for finding seed hits: minimizers (default), maximal unique matches (MUMs) and maximal exact matches (MEMs). Only matches entirely within a node are found. Minimizers (default) are faster and MUM/MEMs can be more sensitive. MUM/MEM modes use [MUMmer4](https://github.com/mummer4/mummer) to find matches between the read and nodes. Use the parameter `--seeds-mum-count n` to use the `n` longest MUMs as seeds (or -1 for all MUMs), and `--seeds-mem-count n` for the `n` longest MEMs (or -1 for all MEMs). Use `--seeds-mxm-length n` to only use matches at least `n` characters long. If you are aligning multiple files to the same graph, use `--seeds-mxm-cache-prefix file_name_prefix` to store the MUM/MEM index to disk for reuse instead of rebuilding it each time.

Alternatively you can use any method to find seed hits and then import the seeds in [.gam format](https://github.com/vgteam/libvgio/blob/master/deps/vg.proto) with the parameter `-s seedfile.gam`. The seeds must be passed as an alignment message, with `path.mapping[0].position` describing the position in the graph, `name` the name of the read and `query_position` the position in the forward strand of the read. Match length (`path.mapping[0].edit[0].from_length`) is only used to order the seeds, with longer matches tried before shorter matches.

Alternatively you can use the parameter `--seeds-first-full-rows` to use the dynamic programming alignment algorithm on the entire first row instead of using seeded alignment. This is very slow except on tiny graphs, and not recommended.

#### Extension

The aligner uses a bitvector banded DP alignment algorithm to extend the seed hits. The DP matrix is calculated inside a certain area (the band), which depends on the extension parameters. Note that "bandwidth" in graph alignment does NOT directly correspond to bandwidth in linear alignment. The bandwidth parameter describes the maximum allowed score difference between the minimum score in a row and a cell, with cells whose score is higher than that falling outside the band. Bandwidth higher than 35 is not recommended for complex graphs due to huge increases in runtime but might work for variation graphs.

The algorithm starts using the initial bandwidth. Should it detect that the alignment is incorrect, it will rewind and rerun with the ramp bandwidth parameter, aligning high-error parts of the read without slowing down alignment in low-error parts. The tangle effort parameter determines how much time the aligner spends inside complex cyclic subgraphs. If the size of the band grows beyond the tangle effort parameter, the aligner will use the current best alignment for the aligned prefix and move forward along the read. This might miss the optimal alignment.

### Parameters

- `-g` input graph. Format .gfa / .vg
- `-f` input reads. Format .fasta / .fastq / .fasta.gz / .fastq.gz. You can input multiple files with `-f file1 -f file2 ...` or `-f file1 file2 ...`
- `-t` number of aligner threads. The program also uses two IO threads in addition to these.
- `-a` output file name. Format .gam or .json
- `--try-all-seeds` extend from all seeds. Normally a seed is not extended if it looks like a false positive.
- `--all-alignments` output all alignments. Normally only a set of non-overlapping partial alignments is returned. Use this to also include partial alignments which overlap each others. This also forces `--try-all-seeds`.
- `--global-alignment` force the read to be aligned end-to-end. Normally the alignment is stopped if the score gets too poor. This forces the alignment to continue to the end of the read regardless of score. If you use this you should do some other filtering on the alignments to remove false alignments.
- `-x` parameter preset. Use `-x vg` for aligning to variation graphs and other simple graphs, and `-x dbg` for aligning to de Bruijn graphs.

Seeding:

- `-s` External seeds. Load seeds from a .gam file. You can input multiple files with `-s file1 -s file2 ...` or `-s file1 file2 ...`
- `--seeds-minimizer-density` For a read of length `n`, use the `arg * n` most unique seeds
- `--seeds-minimizer-length` k-mer size for minimizer seeds
- `--seeds-minimizer-windowsize` Window size for minimizer seeds
- `--seeds-mum-count` MUM seeds. Use the n longest maximal unique matches. -1 for all MUMs
- `--seeds-mem-count` MEM seeds. Use the n longest maximal exact matches. -1 for all MEMs
- `--seeds-mxm-length` MUM/MEM minimum length. Don't use MUMs/MEMs shorter than n
- `--seeds-mxm-cache-prefix` MUM/MEM file cache prefix. Store the MUM/MEM index into disk for reuse. Recommended unless you are sure you won't align to the same graph multiple times
- `--seeds-first-full-rows` Don't use seeds. Instead use the DP alignment on the first row. The runtime depends on the size of the graph so this is very slow. Not recommended

Extension:

- `-b` alignment bandwidth. Unlike in linear alignment, this is the score difference between the minimum score in a row and the score where a cell falls out of the band. Values recommended to be between 1-35.
- `-B` ramp bandwidth. If a read cannot be aligned with the alignment bandwidth, switch to the ramp bandwidth at the problematic location. Values recommended to be between 1-35.
- `-C` tangle effort. Determines how much effort the aligner spends on tangled areas. Higher values use more CPU and memory and have a higher chance of aligning through tangles. Lower values are faster but might return an inoptimal or a partial alignment. Use for complex graphs (eg. de Bruijn graphs of mammalian genomes) to limit the runtime in difficult areas. Values recommended to be between 1'000 - 500'000.
- `--high-memory` high memory mode. Runs a bit faster but uses a LOT more memory
