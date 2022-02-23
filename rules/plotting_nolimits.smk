## We're not using this plot. It is the R version of 'at least one significant.'
# rule plot_sig_enrichment_counts:
#     input:
#         "data/results/cameraPR_nolimits/aggregation/dataId_isSig.tsv"
#     output:
#         "figures/cameraPR_nolimits/at_least_one_significant_facet_R.svg"
#     log:
#         "logs/plot_sig_enrichment_counts.log"
#     conda:
#         "envs/r363.yml"
#     shell:
#         "Rscript scripts/plot_sig_enrichment_counts.R {input} {output} &> {log}"


rule plot_metrics_nolimits:
    input:
        dataId_isSig_file = "data/results/cameraPR_nolimits/aggregation/dataId_isSig.tsv",
        sigTerm_file = "data/results/cameraPR_nolimits/aggregation/sigTermDf.tsv"
    output:
        "figures/cameraPR_nolimits/at_least_one_significant_facetGrid.svg",
        "figures/cameraPR_nolimits/nr_sig_terms_per_user_input.svg",
    params:
        output_dir = lambda wildcards, output: output[0][0:25]
    log:
        "logs/plot_metrics.log"
    conda:
        "envs/py38_plotting.yml"
    shell:
        "python scripts/plot_metrics.py {input.sigTerm_file} {input.dataId_isSig_file} {params.output_dir} &> {log}"



rule plot_metrics_reduced_species_nolimits:
    input:
        dataId_isSig_file = "data/results/cameraPR_nolimits/aggregation/dataId_isSig.tsv",
        sigTerm_file = "data/results/cameraPR_nolimits/aggregation/sigTermDf.tsv"
    output:
        "figures/cameraPR_nolimits/at_least_one_significant_facetGrid_reduced_species.svg",
        "figures/cameraPR_nolimits/nr_sig_terms_per_user_input_reduced_species.svg",
    params:
        output_dir = lambda wildcards, output: output[0][0:25]
    log:
        "logs/plot_metrics_reduced_species.log"
    conda:
        "envs/py38_plotting.yml"
    shell:
        "python scripts/plot_metrics_reduced_species.py {input.sigTerm_file} {input.dataId_isSig_file} {params.output_dir} &> {log}"
