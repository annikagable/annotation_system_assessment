rule downsample_PubMed:
    input:
        effect_file = "data/results/cameraPR/overlap_{_min}-{_max}/effect_size/{dataId}.{taxId}.{db}.tsv"
    output:
        sampled_effect_file = "data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/effect_size/{dataId}.{taxId}.{db}.tsv"
    params:
        sampling_factor = 100
    log:
        "logs/downsample_PubMed/overlap_{_min}-{_max}/{dataId}.{taxId}.{db}.log"
    shell:
        "scripts/downsample_single_PubMed_enrichment.sh {input.effect_file} {params.sampling_factor} {output.sampled_effect_file} &> {log}"


def collect_downsampled_PubMed(wildcards):

    filtered_dedup_dir = checkpoints.apply_deduplication.get().output[1]
    DATAIDS, TAXIDS = glob_wildcards(os.path.join(filtered_dedup_dir,"{dataId}.{taxId}.input.tsv"))

    result_files =  expand(
                        expand("data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/effect_size/{dataId}.{taxId}.{db}.tsv",
                        zip, dataId = DATAIDS, taxId = TAXIDS, allow_missing = True),
                    db = DATABASES, allow_missing = True)
                
    return result_files


rule collect_cameraPR_downsampled_PubMed:
    input:
        collect_downsampled_PubMed
    output:
        "data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/aggregation/aggregated.txt"
    log:
        "logs/collect_cameraPR_downsampled_PubMed.overlap_{_min}-{_max}.log"
    run:
        result_files_string = '\n'.join(input)+'\n'
        with open (output[0], 'w') as f:
            f.write(result_files_string)
            
            
            
rule write_cameraPR_termDf_downsampled_PubMed:
    input:
        enrichment_files_file = "data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/aggregation/aggregated.txt",
        species_taxIds_file = "data/raw/species.v11.0.txt"
    output:
        "data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/aggregation/sigTermDf_alpha1.0.tsv",
        "data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/aggregation/dataId_isSig_alpha1.0.tsv"
    params:
        alpha = 1.0, #get also insignificant terms for plotting term size
        n_grouped_species = 8,
        output_dir = lambda wildcards, output: os.path.dirname(output[0]),
        enrichment_method = "cameraPR"
    log:
        "logs/write_cameraPR_termDf_downsampled_PubMed.overlap_{_min}-{_max}.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/read_enrichment_results.py {input.enrichment_files_file} {input.species_taxIds_file} {params.alpha} {params.n_grouped_species} {params.enrichment_method} {params.output_dir} &> {log}"
        
        
        
rule plot_term_size_v_effect:
    input:
        "data/raw/global_enrichment_annotations/9606.terms_members.tsv",
        "data/results/cameraPR/overlap_{_min}-{_max}/downsampled_PubMed/aggregation/sigTermDf_alpha1.0.tsv"
    output:
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_pval.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_pval_by_database.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_effect_size.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_effect_size_all_terms.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_effect_size_by_database.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_effect_size_FDR.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.term_size_v_effect_size_FDR_by_database.svg",
        "figures/cameraPR/overlap_{_min}-{_max}/term_size/9606.sig_v_insig_term_size.svg",
        "data/results/cameraPR/overlap_{_min}-{_max}/term_size_v_effect_correlations.tsv" 
    log:
        "logs/plot_term_size_v_effect.overlap_{_min}-{_max}.log"
    conda:
        "../envs/py38_plotting.yml"
    shell:
        "python scripts/plot_term_size_v_effect.py {input} {output} >& {log}"

