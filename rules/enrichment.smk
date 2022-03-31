rule ascertain_annotation_is_sorted_by_db:
    input:
        members_file = "data/raw/global_enrichment_annotations/{taxId}.terms_members.tsv"
    output:
        flag_file = "data/interim/annotation_terms/{taxId}.is.sorted"
    shell:
        """
        uniq_etype_count=$(cut -f2 {input} | sort | uniq | wc -l)
        etype_count=$(cut -f2 {input} | uniq | wc -l) 
     
        if [ "$uniq_etype_count" == "$etype_count" ]; then
            touch {output}
        else
            echo "{input} is not sorted by database (aka etype)!"
            exit 1
        fi
        """

rule parse_term_lists:
    input:
        members_file = "data/raw/global_enrichment_annotations/{taxId}.terms_members.tsv",
        flag_file = "data/interim/annotation_terms/{taxId}.is.sorted"
    output:
        expand("data/interim/annotation_terms/{taxId}.term_list.{db}.rds", db = DATABASES, allow_missing = True)
    params:
        output_dir = lambda wildcards, output: output[0][:29]
    log:
        log_file = "logs/parse_term_lists/{taxId}.log"
    conda:
        "../envs/r363.yml"
    script:
        "../scripts/parse_term_lists_by_db.R"


rule run_cameraPR:
    input:
        dedup_id_file = "data/interim/deduplicated_dataIds/{taxId}.tsv",
        infile = "data/interim/filtered_deduplicated_user_inputs/{dataId}.{taxId}.input.tsv", 
        database_file = "data/interim/annotation_terms/{taxId}.term_list.{db}.rds"
    output:
        output_file = "data/results/cameraPR/enrichment/{dataId}.{taxId}.{db}.tsv"
    params:
        min_overlap = 3,
        max_overlap = 200
    log:
        log_file = "logs/run_cameraPR/{dataId}.{taxId}.{db}.log"
    conda:
        "../envs/r363.yml"
    script:
        "../scripts/run_cameraPR.R"
        


# One could potentially do the effect size calculation only once per dataId, i.e. for 
# all databases at the same time. This would be more efficient and faster, but also make
# the input/output handling more complicated. 
# The expand functionality could be used to make a list of input files:
# expand("{dataset}/a.txt", dataset=DATASETS)
    
        
rule get_effect_size:
    input:
        enrichment_file = "data/results/cameraPR/enrichment/{dataId}.{taxId}.{db}.tsv",
        user_input_file = "data/interim/filtered_deduplicated_user_inputs/{dataId}.{taxId}.input.tsv"
    output:
        effect_file = "data/results/cameraPR/effect_size/{dataId}.{taxId}.{db}.tsv"
    log:
        log_file = "logs/get_effect_size/{dataId}.{taxId}.{db}.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/get_effect_size.py {input.enrichment_file} {input.user_input_file} {output} &> {log}"
        

# Collecting cameraPR output only after effect size calculation
def collect_cameraPR_output(wildcards):

    filtered_dedup_dir = checkpoints.apply_deduplication.get().output[1]
    DATAIDS, TAXIDS = glob_wildcards(os.path.join(filtered_dedup_dir,"{dataId}.{taxId}.input.tsv"))

    result_files =  expand(
                expand("data/results/cameraPR/effect_size/{dataId}.{taxId}.{db}.tsv",
                        zip, dataId = DATAIDS, taxId = TAXIDS, allow_missing = True), 
                db = DATABASES)
                
    #result_files = [f for f in result_files if '83332' in f] + [f for f in result_files if 'zScfigOSQ47O' in f]
    #result_files = [f for f in result_files if 'PubMed' not in f]
    return result_files



rule collect_cameraPR:
    input:
        collect_cameraPR_output
    output:
        "data/results/cameraPR/aggregation/aggregated.txt"
    log:
        "logs/collect_cameraPR.log"
    run:
        result_files_string = '\n'.join(input)+'\n'
        with open (output[0], 'w') as f:
            f.write(result_files_string)
            


rule write_cameraPR_termDf:
    input:
        enrichment_files_file = "data/results/cameraPR/aggregation/aggregated.txt",
        species_taxIds_file = "data/raw/species.v11.0.txt"
    output:
        "data/results/cameraPR/aggregation/sigTermDf.tsv",
        "data/results/cameraPR/aggregation/dataId_isSig.tsv"
    params:
        alpha = 0.05,
        n_grouped_species = 8,
        enrichment_method = "cameraPR",
        output_dir = "data/results/cameraPR/aggregation"
    log:
        "logs/write_cameraPR_termDf.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/read_enrichment_results.py {input.enrichment_files_file} {input.species_taxIds_file} {params.alpha} {params.n_grouped_species} {params.enrichment_method} {params.output_dir} &> {log}"





