rule calc_uniq_enriched_genes_per_term_size:
    input:
        sigTermDf_file = "data/results/cameraPR/aggregation/sigTermDf.tsv" 
        hsa_members_file = "data/raw/global_enrichment_annotations/9606.terms_members.tsv"
    params:
        taxId = 9606
    output:
        # union of enriched genes at a specific term size.
        gene_unions_by_term_size_file = "data/results/cameraPR/unique_enriched_proteins/gene_unions_by_term_size.tsv",
        # computed on all possible term sizes => for plotting summary curves
        summary_cum_gene_unions_file = "data/results/cameraPR/unique_enriched_proteins/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv",
#         # computed only on existing term sizes => for plotting individual users curves
#         cum_gene_unions_file= "data/results/cameraPR/unique_enriched_proteins/cumulative_gene_unions_by_term_size.tsv"
    log:
        "logs/calc_uniq_enriched_genes_per_term_size.log"



