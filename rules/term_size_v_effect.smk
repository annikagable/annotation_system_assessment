#FNAMES = glob_wildcards("data/results/cameraPR/effect_size/{fname}").fname

rule downsample_PubMed:
    input:
        "data/results/cameraPR/aggregation/aggregated.txt"
    output:
        expand("data/results/cameraPR/effect_size_downsampled_PubMed/{fname}", fname = glob_wildcards("data/results/cameraPR/effect_size/{fname}").fname)
    params:
        input_dir = "data/results/cameraPR/effect_size",
        sampling_factor = 100,
        output_dir = "data/results/cameraPR/effect_size_downsampled_PubMed"
    log:
        "logs/downsample_PubMed.log"
    shell:
        "scripts/downsample_PubMed_enrichments.sh {params.input_dir} {params.sampling_factor} {params.output_dir} &> {log}"
        
        
        
rule collect_cameraPR_downsampled_PubMed:
    input:
        expand("data/results/cameraPR/effect_size_downsampled_PubMed/{fname}",
               fname = glob_wildcards("data/results/cameraPR/effect_size/{fname}").fname)
    output:
        "data/results/cameraPR/aggregation_downsampled_PubMed/aggregated.txt"
    log:
        "logs/collect_cameraPR_downsampled_PubMed.log"
    run:
        result_files_string = '\n'.join(input)+'\n'
        with open (output[0], 'w') as f:
            f.write(result_files_string)
            
            
            
rule write_cameraPR_termDf_downsampled_PubMed:
    input:
        enrichment_files_file = "data/results/cameraPR/aggregation_downsampled_PubMed/aggregated.txt",
        species_taxIds_file = "data/raw/species.v11.0.txt"
    output:
        "data/results/cameraPR/aggregation_downsampled_PubMed/sigTermDf.tsv",
        "data/results/cameraPR/aggregation_downsampled_PubMed/dataId_isSig.tsv"
    params:
        alpha = 1,
        n_grouped_species = 8,
        output_dir = "data/results/cameraPR/aggregation_downsampled_PubMed",
        enrichment_method = "cameraPR"
    log:
        "logs/write_cameraPR_termDf_downsampled_PubMed.log"
    conda:
        "../envs/py38_sklearn.yml"
    shell:
        "python scripts/read_enrichment_results.py {input.enrichment_files_file} {input.species_taxIds_file} {params.alpha} {params.n_grouped_species} {params.enrichment_method} {params.output_dir} &> {log}"
        
        
        
rule plot_term_size_v_effect:
    input:
        "data/raw/global_enrichment_annotations/9606.terms_members.tsv",
        "data/results/cameraPR/aggregation_downsampled_PubMed/sigTermDf.tsv"
    output:
        "figures/cameraPR/9606.term_size_v_pval.svg",
        "figures/cameraPR/9606.term_size_v_pval_by_database.svg",
        "figures/cameraPR/9606.term_size_v_effect_size.svg",
        "figures/cameraPR/9606.term_size_v_effect_size_all_terms.svg",
        "figures/cameraPR/9606.term_size_v_effect_size_by_database.svg",
        "figures/cameraPR/9606.term_size_v_effect_size_FDR.svg",
        "figures/cameraPR/9606.term_size_v_effect_size_FDR_by_database.svg",
        "figures/cameraPR/9606.sig_v_insig_term_size.svg",
        "data/results/cameraPR/term_size_v_effect_correlations.tsv" 
    log:
        "logs/plot_term_size_v_effect.log"
    conda:
        "../envs/py38_plotting.yml"
    shell:
        "python scripts/plot_term_size_v_effect.py {input} {output} >& {log}"
