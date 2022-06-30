rule filter_inputs:
    input:
        config["user_input_dir"]
    output:
        "data/interim/filtered_user_inputs.tsv",
    conda:
        "../envs/py38_sklearn.yml"
    params:
        search_for = config["input_file_pattern"],
        min_finite_values = config["min_finite_values"],
        min_unique_values = config["min_unique_values"],
        max_ratio_top_value_to_input_size = config["max_ratio_top_value_to_input_size"]
    log:
        "logs/filter_inputs.log"
    shell:
        "python scripts/filter_user_inputs.py --search_for {params.search_for} --min_finite_values {params.min_finite_values} --min_unique_values {params.min_unique_values} --max_ratio_top_value_to_input_size {params.max_ratio_top_value_to_input_size} --out_file {output} {input} &> {log}"


checkpoint separate_inputs_by_species:
    input:
        "data/interim/filtered_user_inputs.tsv"
    output:
        directory("data/interim/species_matrices") # this is where e.g. 9606.tsv files will go for each species
    conda:
        "../envs/py38_sklearn.yml"
    log:
        "logs/separate_inputs_by_species.log"
    shell:
        "python scripts/separate_inputs_by_species.py {input} {output} &> {log}"


rule deduplicate_similar_inputs:
    input:
        species_matrix_dir = "data/interim/species_matrices",
        species_matrix_file = "data/interim/species_matrices/{taxId}.tsv"
    params:
        parallelization_threshold = config["parallelization_threshold"],
        min_overlap = config["min_overlap"],
        r2_threshold = config["r2_threshold"],
        symm_diff_threshold = config["symm_diff_threshold"],
        metric_matrix_dir = lambda wildcards, output: os.path.dirname(output[1]),
    output:
        dedup_id_file     = "data/interim/deduplicated_dataIds/{taxId}.tsv",
        spearman_file     = "data/interim/metric_matrices/{taxId}/spearman.overlap_"+config["min_overlap"]+".tsv",
        spearman_abs_file = "data/interim/metric_matrices/{taxId}/spearman_abs.overlap_"+config["min_overlap"]+".tsv",
        sym_diff_file     = "data/interim/metric_matrices/{taxId}/sym_diff.tsv",
        sl_cluster_file   = "data/interim/metric_matrices/{taxId}/single_linkage_clusters.tsv"
    threads:
        workflow.cores
    conda:
        "../envs/py38_sklearn.yml"
    log:
        "logs/deduplication/{taxId}.log"
    shell:
        "python scripts/deduplicate_similar_inputs.py {input.species_matrix_file} {threads} {params} {output.dedup_id_file}  &> {log}"




def expand_deduplicated_dataId_files(wildcards):

    species_matrix_dir = checkpoints.separate_inputs_by_species.get().output[0]
    TAXIDS = glob_wildcards(os.path.join(species_matrix_dir, "{taxId}.tsv")).taxId
    dedup_dataId_files = expand("data/interim/deduplicated_dataIds/{taxId}.tsv", taxId = TAXIDS)
    return dedup_dataId_files

checkpoint apply_deduplication:
    input:
        filtered_user_inputs_file = "data/interim/filtered_user_inputs.tsv",
        deduplicate_dataId_files =  expand_deduplicated_dataId_files
    output:
        table_out_file = "data/interim/filtered_deduplicated_user_inputs.tsv",
        out_dir = directory("data/interim/filtered_deduplicated_user_inputs"),
    conda:
        "../envs/py38_sklearn.yml"
    log:
        "logs/deduplication/filter_deduplicated.log"
    script:
        "../scripts/filter_deduplicated.py"


rule report_filtered_and_deduplicated_count:
    input:
        table_out_file = "data/interim/filtered_deduplicated_user_inputs.tsv"
    output:
        report_file = "data/results/filtering_report.txt",
        by_taxId_file = "data/results/deduplicated_user_submission_counts_by_taxId.tsv"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/report_number_of_filtered_deduplicated_inputs.py {input.table_out_file} {output.by_taxId_file} &> {output.report_file}"



