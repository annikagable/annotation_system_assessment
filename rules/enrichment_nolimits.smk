
rule run_cameraPR_nolimits:
    input:
        dedup_id_file = "data/interim/deduplicated_dataIds/{taxId}.tsv",
        infile = "data/interim/filtered_deduplicated_user_inputs/{dataId}.{taxId}.input.tsv", 
        database_file = "data/interim/annotation_terms/{taxId}.term_list.{db}.rds"
    output:
        output_file = "data/results/cameraPR_nolimits/enrichment/{dataId}.{taxId}.{db}.tsv"
    log:
        log_file = "logs/run_cameraPR_nolimits/{dataId}.{taxId}.{db}.log"
    conda:
        "../envs/r363.yml"
    script:
        "../scripts/run_cameraPR_nolimits.R"

        
rule get_effect_size_nolimits:
    input:
        enrichment_file = "data/results/cameraPR_nolimits/enrichment/{dataId}.{taxId}.{db}.tsv",
        user_input_file = "data/interim/filtered_deduplicated_user_inputs/{dataId}.{taxId}.input.tsv"
    output:
        effect_file = "data/results/cameraPR_nolimits/effect_size/{dataId}.{taxId}.{db}.tsv"
    log:
        log_file = "logs/get_effect_size_nolimits/{dataId}.{taxId}.{db}.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/get_effect_size.py {input.enrichment_file} {input.user_input_file} {output} &> {log}"
        

def collect_cameraPR_output_nolimits(wildcards):

    filtered_dedup_dir = checkpoints.apply_deduplication.get().output[1]
    DATAIDS, TAXIDS = glob_wildcards(os.path.join(filtered_dedup_dir,"{dataId}.{taxId}.input.tsv"))

    result_files =  expand(
                expand("data/results/cameraPR_nolimits/effect_size/{dataId}.{taxId}.{db}.tsv",
                        zip, dataId = DATAIDS, taxId = TAXIDS, allow_missing = True), 
                db = DATABASES)
                
    return result_files



rule collect_cameraPR_nolimits:
    input:
        collect_cameraPR_output_nolimits
    output:
        "data/results/cameraPR_nolimits/aggregation/aggregated.txt"
    log:
        "logs/collect_cameraPR_nolimits.log"
    run:
        result_files_string = '\n'.join(input)+'\n'
        with open (output[0], 'w') as f:
            f.write(result_files_string)
            


rule write_cameraPR_termDf_nolimits:
    input:
        enrichment_files_file = "data/results/cameraPR_nolimits/aggregation/aggregated.txt",
        species_taxIds_file = "data/raw/species.v11.0.txt"
    output:
        "data/results/cameraPR_nolimits/aggregation/sigTermDf.tsv",
        "data/results/cameraPR_nolimits/aggregation/dataId_isSig.tsv"
    params:
        alpha = 0.05,
        n_grouped_species = 8,
        enrichment_method = "cameraPR",
        output_dir = lambda wildcards, output: output[0][0:42],
    log:
        "logs/write_cameraPR_termDf_nolimits.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/read_enrichment_results.py {input.enrichment_files_file} {input.species_taxIds_file} {params.alpha} {params.n_grouped_species} {params.enrichment_method} {params.output_dir} &> {log}"





