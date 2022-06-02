import os
import glob

#configfile: "config.yaml"
SPECIES_SUBSETS = ["all", "reduced"]
OVERLAP_THRESHOLDS = [(3,"Inf"), (0, "Inf")]
MINS, MAXS = zip(*OVERLAP_THRESHOLDS)

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
USER_INPUT_DIR = "data/raw/string_11_0_gene_value_consolidated_input"  #"data/raw/example_user_queries" #not used yet
#ANNOTATION_DIR = "data/raw/global_enrichment_annotations" #not used yet


# significance threshold for enrichment 
ALPHA = 0.05

wildcard_constraints:
    _min="[a-zA-Z0-9]+",
    _max="[a-zA-Z0-9]+",
    lo = "[a-zA-Z0-9]+",
    hi = "[a-zA-Z0-9]+"

#DATAIDS, TAXIDS = [["zDPnSS4z4TPV", "z4WZvoMHP74T"], ["9606", "10090"]]

def expand_plots_by_overlap_thresholds(wildcards):
    
    MINS, MAXS = zip(*OVERLAP_THRESHOLDS)
    
    enr_plots = [
            "at_least_one_significant_facetGrid", 
            "nr_sig_terms_per_user_input",
            ]

    term_size_plots = [
        "9606.term_size_v_pval",
        "9606.term_size_v_pval_by_database",
        "9606.term_size_v_effect_size",
        "9606.term_size_v_effect_size_all_terms",
        "9606.term_size_v_effect_size_by_database",
        "9606.term_size_v_effect_size_FDR",
        "9606.term_size_v_effect_size_FDR_by_database",
        "9606.sig_v_insig_term_size"
        ]

    enrichment_plots = expand(
                           expand("figures/cameraPR/overlap_{_min}-{_max}/enriched_terms/{species_subset}_species/{plot}.svg", 
                           zip, _min = MINS, _max = MAXS, allow_missing = True),
                       species_subset = SPECIES_SUBSETS, plot = enr_plots)

    term_size_v_pval_plots = expand(
                                expand("figures/cameraPR/overlap_{_min}-{_max}/term_size/{term_size_plot}.svg",
                                zip, _min = MINS, _max = MAXS, allow_missing = True),
                             term_size_plot = term_size_plots)

    redundancy_and_novelty_plots = expand(
                                        expand("figures/cameraPR/overlap_{_min}-{_max}/redundancy_and_novelty/{species_subset}_species/{plot}",
                                        zip, _min = MINS, _max = MAXS, allow_missing = True),
                                   species_subset = SPECIES_SUBSETS, plot = ["nr_nonredundant_sig_terms_per_user_input.svg",
                                                                             "at_least_one_novel_significant.svg"])

    unique_enriched_genes_plots = expand("figures/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/9606/enrichment_specificity.svg",
            #"figures/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/9606/enrichment_specificity_legend.svg",
                                  zip, _min=MINS, _max = MAXS)

    result_files = enrichment_plots + term_size_v_pval_plots + redundancy_and_novelty_plots + unique_enriched_genes_plots
    
    return result_files



rule all:
    input:
        expand_plots_by_overlap_thresholds,
        "figures/input_analysis/input_count_and_input_size_by_species_group.svg",

        #"figures/database_stats/9606.database_stats.term_size_11708-Inf.pdf",
        "figures/database_stats/9606.database_stats.term_size_0-250.pdf",
        "figures/database_stats/9606.database_stats.term_size_0-Inf.pdf",


        "data/results/filtering_report.txt",
        "data/results/deduplicated_user_submission_counts_by_taxId.tsv",

        expand("data/results/cameraPR/overlap_{_min}-{_max}/term_size_v_effect_correlations.tsv",
                zip, _min = MINS, _max = MAXS) 


include:
    "rules/filter_and_deduplicate.smk"
include:
    "rules/enrichment.smk"
include:
    "rules/term_size_v_effect.smk"
include:
    "rules/user_input_stats.smk"
include:
    "rules/redundancy_and_novelty.smk"
include:
    "rules/enrichment_specificity.smk"
include:
    "rules/database_stats.smk"
