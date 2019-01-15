A Snakemake pipeline for error correcting long reads based on short reads.

Installation:

- Install snakemake https://snakemake.readthedocs.io/en/stable/
- Install lighter https://github.com/mourisl/Lighter
- Install bcalm2 https://github.com/GATB/bcalm
- Install the aligner and `make all`

Running:

- Edit the parameters in `config.yaml`
- `snakemake all`
- The corrected reads will be in the output folder in `corrected.fasta`
