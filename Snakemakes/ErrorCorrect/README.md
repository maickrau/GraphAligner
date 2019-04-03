A Snakemake pipeline for error correcting long reads based on short reads.

Installation:

- Install snakemake https://snakemake.readthedocs.io/en/stable/
- Install lighter https://github.com/mourisl/Lighter
- Install bcalm2 https://github.com/GATB/bcalm
- Build the aligner with `make all`. Note that the pipeline uses scripts which are not included in the bioconda release so you have to compile it from source.

Running:

- Edit the parameters in `config.yaml`
- `snakemake all`
- The corrected reads will be in the output folder:
  - `corrected.fasta` has the reads with the aligned sequence replaced by the alignment path. Uppercase sequence are corrected and lowercase are uncorrected.
  - `corrected_clipped.fasta` has the reads cut across non-corrected parts. All sequence is corrected and the read name contains the position in the original read.
