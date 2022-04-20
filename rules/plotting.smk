rule plot_metrics:
    input:
        dataId_isSig_file = "data/results/cameraPR/overlap_{_min}-{_max}/aggregation/dataId_isSig.tsv",
        sigTerm_file = "data/results/cameraPR/overlap_{_min}-{_max}/aggregation/sigTermDf.tsv"
    output:
        "figures/cameraPR/overlap_{_min}-{_max}/at_least_one_significant_facetGrid.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/nr_sig_terms_per_user_input.svg"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output[0])
    log:
        "logs/plot_metrics/overlap_{_min}-{_max}.log"
    conda:
        "../envs/py38_plotting.yml"
    shell:
        "python scripts/plot_metrics.py {input.sigTerm_file} {input.dataId_isSig_file} {params.output_dir} &> {log}"



rule plot_metrics_reduced_species:
    input:
        dataId_isSig_file = "data/results/cameraPR/aggregation/dataId_isSig.tsv",
        sigTerm_file = "data/results/cameraPR/aggregation/sigTermDf.tsv"
    output:
        "figures/cameraPR/at_least_one_significant_facetGrid_reduced_species.svg",
        "figures/cameraPR/nr_sig_terms_per_user_input_reduced_species.svg"
    params:
        output_dir = lambda wildcards, output: output[0][0:16]
    log:
        "logs/plot_metrics_reduced_species.log"
    conda:
        "../envs/py38_plotting.yml"
    shell:
        "python scripts/plot_metrics_reduced_species.py {input.sigTerm_file} {input.dataId_isSig_file} {params.output_dir} &> {log}"
