rule calc_uniq_enriched_genes_per_term_size:
    input:
        sigTermDf = "data/results/cameraPR/overlap_3-200/aggregation/sigTermDf_alpha"+str(ALPHA)+".tsv", 
        species_members = "data/raw/global_enrichment_annotations/9606.terms_members.tsv"
    params:
        taxId = lambda wildcards,  input: '.'.split(os.path.basename(input[1]))[0] #9606
    output:
        # union of enriched genes at a specific term size.
        gene_unions_by_term_size = "data/results/cameraPR/overlap_3-200/unique_enriched_genes/gene_unions_by_term_size.tsv",
        # computed on all possible term sizes => for plotting summary curves
        summary_cum_gene_unions = "data/results/cameraPR/overlap_3-200/unique_enriched_genes/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv",
        # computed only on existing term sizes => for plotting individual users curves
        #cum_gene_unions = "data/results/cameraPR/overlap_3-200/unique_enriched_genes/cumulative_gene_unions_by_term_size.tsv"
    log:
        "logs/calc_uniq_enriched_genes_per_term_size.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/calc_uniq_enriched_genes_per_term_size.py {input.sigTermDf} {input.species_members} {params.taxId} {output.gene_unions_by_term_size} {output.summary_cum_gene_unions} &> {log}"



