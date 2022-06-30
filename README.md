# Filtering and deduplication of user submissions to STRING global enrichment

This is a part of the snakemake pipeline annotation_system_assessment, as of June 29th, 2022.

The filtering and duplication settings here are what was being used for the manuscript
"Systematic assessment of pathway databases, based on a diverse collection of user-submitted experiments".

To use:

Create a conda environment for running the pipeline:
`conda env create -f envs/snakemake.yml`

Activate the environment:
`conda activate snakemake`

Run the pipeline:
`snakemake --use-conda --cores 10`

In order to change the input directory and other parameters, edit the `config.yaml`.

