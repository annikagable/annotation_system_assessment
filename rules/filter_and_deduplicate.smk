#DATAIDS = ["zz26G5go7Urs"]
#DATABASES = ["KEGG", "STRINGclusters"]

# import os
# import glob




rule filter_inputs:
    input:
        USER_INPUT_FOLDER
    output:
        "data/interim/filtered_user_inputs.tsv",
    conda:
        "../envs/py38_sklearn.yml"
    params:
        search_for = '*.input.normal.txt',
        min_finite_values = '500',
        min_unique_values = '10',
        max_ratio_top_value_to_input_size = '0.8'
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
    shell:
        "python scripts/separate_inputs_by_species.py {input} data/interim/species_matrices"


rule deduplicate_similar_inputs:
    input:
        species_matrix_dir = "data/interim/species_matrices",
        species_matrix_file = "data/interim/species_matrices/{taxId}.tsv"
    params:
        metric_matrix_dir = lambda wildcards, output: os.path.dirname(output[1]),
        min_overlap = lambda wildcards, output: os.path.basename(output[1]).split('_')[1].split('.')[0]
    output:
        dedup_id_file     = "data/interim/deduplicated_dataIds/{taxId}.tsv",
        spearman_file     = "data/interim/metric_matrices/{taxId}/spearman.overlap_100.tsv",
        spearman_abs_file = "data/interim/metric_matrices/{taxId}/spearman_abs.overlap_100.tsv",
        sym_diff_file     = "data/interim/metric_matrices/{taxId}/sym_diff.tsv",
        sl_cluster_file   = "data/interim/metric_matrices/{taxId}/single_linkage_clusters.tsv"
    threads:
        workflow.cores
    conda:
        "../envs/py38_sklearn.yml"
    log:
        "logs/deduplication/{taxId}.log"
    shell:
        "python scripts/deduplicate_similar_inputs.py {input.species_matrix_file} {params.metric_matrix_dir} {threads} {params.min_overlap} {output.dedup_id_file}  &> {log}"




def expand_deduplicated_dataId_files(wildcards):

    species_matrix_dir = checkpoints.separate_inputs_by_species.get().output[0]
    TAXIDS = glob_wildcards(os.path.join(species_matrix_dir, "{taxId}.tsv")).taxId
    #print(TAXIDS)
    dedup_dataId_files = expand("data/interim/deduplicated_dataIds/{taxId}.tsv", taxId = TAXIDS)
    return dedup_dataId_files

checkpoint apply_deduplication:
    input:
        filtered_user_inputs_file = "data/interim/filtered_user_inputs.tsv",
        deduplicate_dataId_files =  expand_deduplicated_dataId_files
    output:
        table_out_file = "data/interim/filtered_deduplicated_user_inputs.tsv",
        out_dir = directory("data/interim/filtered_deduplicated_user_inputs"),
        #out_files = "data/interim/filtered_deduplicated_user_inputs/{dataId}.{taxId}.input.tsv"
    conda:
        "../envs/py38_sklearn.yml"
    log:
        "logs/deduplication/filter_deduplicated.log"
    script:
        "../scripts/filter_deduplicated.py"


# def collect_filtered_deduplicated_inputs(wildcards):
# 
#     ## collecting the taxIds of the species_matrices
#     unique_TAXIDS, = glob_wildcards(os.path.join(checkpoints.separate_inputs_by_species.get().output[0], "{taxId}.tsv"))
# 
#     DATAIDS = []
#     TAXIDS = []
#     for _taxId in unique_TAXIDS:
#         output_deduplicate = f"data/interim/deduplicated_dataIds/{_taxId}.tsv"
#         with open(output_deduplicate, 'r') as f:
#             lines = f.readlines()
#         DATAIDS += [l.strip() for l in lines]
#         TAXIDS += [_taxId] * len(lines)
# 
#     checkpoint_output = checkpoints.apply_deduplication.get().output[1]
#     result_files = expand(os.path.join(checkpoint_output,"{dataId}.{taxId}.input.tsv"),
#                           zip, 
#                           dataId = DATAIDS,
#                           taxId = TAXIDS)
#     return result_files
# 
# 
# rule all_filtering:
#     input:
#         #"data/interim/filtered_deduplicated_user_inputs.tsv"
#         #"data/interim/filtered_deduplicated_user_input_files.tsv",
#         collect_filtered_deduplicated_inputs



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



rule enumerate_user_inputs:
    input:
        "data/raw/string_11_0_gene_value_consolidated_input"
    output:
        "data/results/user_input_enumeration.tsv"
    shell:
        """basename -a $(ls {USER_INPUT_FOLDER}/*.input.normal.txt) | cut -d'.' -f1 | awk 'BEGIN{{OFS="\t"}}{{print $0,NR}}' > {output}"""

