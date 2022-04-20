rule plot_metrics:
    input:
        dataId_isSig_file = "data/results/cameraPR/overlap_{_min}-{_max}/aggregation/dataId_isSig.tsv",
        sigTerm_file = "data/results/cameraPR/overlap_{_min}-{_max}/aggregation/sigTermDf.tsv"
    output:
        "figures/cameraPR/overlap_{_min}-{_max}/{species_subset}_species/at_least_one_significant_facetGrid.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/{species_subset}_species/nr_sig_terms_per_user_input.svg"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output[0]),
        species_subset = "{species_subset}"
    log:
        "logs/plot_metrics/overlap_{_min}-{_max}.species_{species_subset}.log"
    conda:
        "../envs/py38_plotting.yml"
    shell:
        "python scripts/plot_metrics.py {input.sigTerm_file} {input.dataId_isSig_file} {params.output_dir} {params.species_subset} &> {log}"
