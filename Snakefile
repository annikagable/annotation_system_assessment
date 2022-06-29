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
        "data/results/filtering_report.txt",
        "data/results/deduplicated_user_submission_counts_by_taxId.tsv",
        


include:
    "rules/filter_and_deduplicate.smk"
