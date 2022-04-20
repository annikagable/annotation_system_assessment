import os
import glob

#configfile: "config.yaml"
SPECIES_SUBSETS = ["all", "reduced"]
OVERLAP_THRESHOLDS = [(3,200), (0, "Inf")]
DATABASES = ['STRINGclusters',
             #'Reactome',
             #'PubMed',
             #'Pfam',
             #'InterPro',
             'SMART']
             #'KEGG',
             #'UniProt',
             #'GO_MF',
             #'GO_CC',
             #'GO_BP']


#DATAIDS, TAXIDS = [["zDPnSS4z4TPV", "z4WZvoMHP74T"], ["9606", "10090"]]

def expand_plots_by_overlap_thresholds(wildcards):
    MINS, MAXS = zip(*OVERLAP_THRESHOLDS)
    #plots = 
    facetGrid_plots = expand(expand("figures/cameraPR/overlap_{_min}-{_max}/{species_subset}_species/at_least_one_significant_facetGrid.svg", zip, _min = MINS, _max = MAXS,allow_missing = True),species_subset = SPECIES_SUBSETS)
    #sig_terms_plots = expand("figures/cameraPR/overlap_{_min}-{_max}/{species_subset}/nr_sig_terms_per_user_input.svg",        zip, _min = MINS, _max = MAXS)
    result_files = facetGrid_plots# + sig_terms_plots
    return result_files

rule all:
    input:
        expand_plots_by_overlap_thresholds,

        "data/results/filtering_report.txt",
        "data/results/deduplicated_user_submission_counts_by_taxId.tsv",

        # Shouldn't need these at all 
        #"data/results/cameraPR/aggregation/aggregated.txt",
        #"data/results/cameraPR/aggregation_downsampled_PubMed/aggregated.txt",
        #"data/results/cameraPR_nolimits/aggregation/aggregated.txt",
        #"data/results/cameraPR_nolimits/aggregation_downsampled_PubMed/aggregated.txt",

        #"figures/cameraPR/at_least_one_significant_facetGrid.svg",
        #"figures/cameraPR/nr_sig_terms_per_user_input.svg",
        #"figures/cameraPR/at_least_one_significant_facetGrid_reduced_species.svg",
        #"figures/cameraPR/nr_sig_terms_per_user_input_reduced_species.svg",

        #"figures/cameraPR/9606.term_size_v_pval.svg",
        #"figures/cameraPR/9606.term_size_v_pval_by_database.svg",
        #"figures/cameraPR/9606.term_size_v_effect_size.svg",
        #"figures/cameraPR/9606.term_size_v_effect_size_by_database.svg",

        #"figures/cameraPR/redundancy_and_novelty/all_species/nr_sig_terms_per_user_input.svg",
        #"figures/cameraPR/redundancy_and_novelty/reduced_species/nr_sig_terms_per_user_input.svg",
        #"figures/cameraPR/redundancy_and_novelty/all_species/at_least_one_significant.svg",
        #"figures/cameraPR/redundancy_and_novelty/reduced_species/at_least_one_significant.svg",

        #"figures/cameraPR_nolimits/at_least_one_significant_facetGrid.svg",
        #"figures/cameraPR_nolimits/nr_sig_terms_per_user_input.svg",
        #"figures/cameraPR_nolimits/at_least_one_significant_facetGrid_reduced_species.svg",
        #"figures/cameraPR_nolimits/nr_sig_terms_per_user_input_reduced_species.svg",

        #"figures/cameraPR_nolimits/9606.term_size_v_pval.svg",
        #"figures/cameraPR_nolimits/9606.term_size_v_pval_by_database.svg",
        #"figures/cameraPR_nolimits/9606.term_size_v_effect_size.svg",
        #"figures/cameraPR_nolimits/9606.term_size_v_effect_size_by_database.svg",

        "figures/input_analysis/input_count_and_input_size_by_species_group.svg",


include:
    "rules/filter_and_deduplicate.smk"
include:
    "rules/enrichment.smk"
include:
    "rules/plotting.smk"
include:
    "rules/term_size_v_effect.smk"

#include:
#    "rules/enrichment_nolimits.smk"
#include:
#    "rules/plotting_nolimits.smk"
#include:
#    "rules/term_size_v_effect_nolimits.smk"

include:
    "rules/plot_user_input_stats.smk"
include:
    "rules/redundancy_and_novelty.smk"

