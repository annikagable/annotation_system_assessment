import os
import glob


DATABASES = ['STRINGclusters',
             'Reactome',
             'PubMed',
             'Pfam',
             'InterPro',
             'SMART',
             'KEGG',
             'UniProt',
             'GO_MF',
             'GO_CC',
             'GO_BP']
#USER_INPUT_DIR = "data/raw/string_11_0_gene_value_consolidated_input"  
#USER_INPUT_DIR = "data/raw/STRING_global_enrichment_example_inputs"
USER_INPUT_DIR = "data/raw/deduplication_test_inputs"

rule all:
    input:
        # This is where the desired output files are listed.

        # Do filtering and deduplication only
        "data/results/filtering_report.txt",
        "data/results/deduplicated_user_submission_counts_by_taxId.tsv",
        
        # Also plot the user input stats
        "figures/input_analysis/input_count_and_input_size_by_species_group.svg",
        "figures/input_analysis/input_count_and_input_size_for_other_species.svg",
        "figures/input_analysis/histogram_human_user_inputs_per_protein.svg"



include:
    "rules/filter_and_deduplicate.smk",
    "rules/user_input_stats.smk"
