rule get_proteome_sizes:
    input:
        "data/raw/protein.info.v11.0.txt.gz"
    output:
        "data/interim/taxId_to_proteome_size.tsv"
    log:
        "logs/get_proteome_sizes.log"
    shell:
        "scripts/get_proteome_sizes.sh {input} > {output} 2> {log}" 


rule plot_user_input_stats:
    input:
        "data/interim/filtered_deduplicated_user_inputs.tsv",
        "data/raw/species.v11.0.txt",
        "data/raw/proteins_to_shorthands.v11.tsv",
        "data/raw/9606.protein.info.v11.0.txt.gz",
        rules.get_proteome_sizes.output
    output:
        "figures/input_analysis/input_count_and_input_size_by_species_group.svg",
        "figures/input_analysis/input_count_and_input_size_for_other_species.svg",
        "figures/input_analysis/histogram_human_user_inputs_per_protein.svg"
    log:
        "logs/plot_user_input_stats.log"
    conda:
        "../envs/py38_plotting.yml"
    shell:
        "python scripts/plot_user_input_statistics.py {input} {output} &> {log}"
        
        
