import os

configfile:
    "config.yaml"

rule all:
    input:
        # This is where the desired output files are listed.
        # Snakemake will create all intermediate files as necessary including, 
        # in this case, the actual filtered deduplicated user inputs
        "data/results/filtering_report.txt",
        "data/results/deduplicated_user_submission_counts_by_taxId.tsv"


include:
    "rules/filter_and_deduplicate.smk"
