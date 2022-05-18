rule get_term_coverage:
    input:
        proteins_to_shorthands_file = "data/raw/proteins_to_shorthands.v11.tsv",
        protein_info_file = "data/raw/{taxId}.protein.info.v11.0.txt.gz",
        terms_members_file = "data/raw/global_enrichment_annotations/{taxId}.terms_members.tsv"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output[0]),
        taxId = "{taxId}",
        max_term_size =  lambda wildcards, output: os.path.basename(output[0]).split(".")[0].split("-")[1]
    output:
        term_coverage_file = "data/results/database_stats/{taxId}/term_coverage_per_gene_0-{max_term_size}.tsv",
        genome_coverage_file =  "data/results/database_stats/{taxId}/genome_coverage_0-{max_term_size}.tsv",
        percent_genome_coverage_file = "data/results/database_stats/{taxId}/percent_genome_coverage_0-{max_term_size}.tsv"
    log:
        "logs/get_term_coverage/{taxId}.0-{max_term_size}.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/get_term_coverage.py {input} {params.output_dir} {params.taxId} {params.max_term_size} &> {log}"


rule plot_database_stats:
    input:
        term_coverage_file = "data/results/database_stats/{taxId}/term_coverage_per_gene_0-{max_term_size}.tsv",
        genome_coverage_file =  "data/results/database_stats/{taxId}/genome_coverage_0-{max_term_size}.tsv",
        percent_genome_coverage_file = "data/results/database_stats/{taxId}/percent_genome_coverage_0-{max_term_size}.tsv",
        species_members_file= "data/raw/global_enrichment_annotations/{taxId}.terms_members.tsv"
        
    output:
        database_stats_file = "figures/database_stats/{taxId}.database_stats_0-{max_term_size}.pdf"
    log:
        "logs/plot_database_stats/{taxId}.0-{max_term_size}.log"
    conda:
        "../envs/py38_mpl35.yml"
    shell:
        "python scripts/plot_database_stats.py {input} {output} &> {log}"
        
        
