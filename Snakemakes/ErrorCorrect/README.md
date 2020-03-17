A Snakemake pipeline for error correcting long reads based on short reads.

Installation:

- Install snakemake, lighter, bcalm2 and GraphAligner: `conda install -c bioconda snakemake lighter bcalm graphaligner`
- Download the bcalm2 GFA conversion script from https://github.com/GATB/bcalm/blob/master/scripts/convertToGFA.py

Running:

- Save `Snakefile` and `config.yaml` to a folder
- Edit the parameters in `config.yaml`
- Run `snakemake --cores 8 all` (you can use more than 8 cores)
- The corrected reads will be in the output folder:
  - `corrected.fasta` has the reads with the aligned sequence replaced by the alignment path. Uppercase sequence are corrected and lowercase are uncorrected.
  - `corrected_clipped.fasta` has the reads cut across non-corrected parts. All sequence is corrected and the read name contains the position in the original read.
