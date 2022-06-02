rule calc_uniq_enriched_genes_per_term_size:
    input:
        sigTermDf_file = "data/results/cameraPR/overlap_{_min}-{_max}/aggregation/sigTermDf_alpha"+str(ALPHA)+".tsv", 
        species_members_file = "data/raw/global_enrichment_annotations/{taxId}.terms_members.tsv"
    params:
        taxId = "{taxId}"
    output:
        # union of enriched genes at a specific term size.
        #gene_unions_by_term_size = "data/results/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/gene_unions_by_term_size.tsv",
        # computed on all possible term sizes => for plotting summary curves
        summary_cum_gene_unions_file = "data/results/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/{taxId}/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv",
        # computed only on existing term sizes => for plotting individual users curves
        cum_gene_unions_file = "data/results/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/{taxId}/cumulative_gene_unions_by_term_size.tsv"
    log:
        "logs/calc_uniq_enriched_genes_per_term_size/{taxId}.overlap_{_min}-{_max}.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/calc_uniq_enriched_genes_per_term_size.py {input.sigTermDf_file} {input.species_members_file} {params.taxId} {output.cum_gene_unions_file} {output.summary_cum_gene_unions_file} &> {log}"

        
rule plot_enrichment_specificity:
    input:
        summary_cum_gene_unions_file = "data/results/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/{taxId}/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv",
        cum_gene_unions_file = "data/results/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/{taxId}/cumulative_gene_unions_by_term_size.tsv"
    output:
        enrichment_specificity_file = "figures/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/{taxId}/enrichment_specificity.pdf",
        enrichment_specificity_legend_file = "figures/cameraPR/overlap_{_min}-{_max}/unique_enriched_genes/{taxId}/enrichment_specificity_legend.pdf"
    log:
        "logs/plot_enrichment_specificity/{taxId}.overlap_{_min}-{_max}.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/plot_enrichment_specificity.py {input.cum_gene_unions_file} {input.summary_cum_gene_unions_file} {output.enrichment_specificity_file} {output.enrichment_specificity_legend_file}"

                 


