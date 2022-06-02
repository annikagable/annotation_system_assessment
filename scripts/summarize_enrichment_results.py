# To be added to snakemake
# Get the median enrichment counts per user input for human, mouse, and other organisms, for the result summary table.

import pandas as pd
import numpy as np
import utils

db_order = utils.dbColors.index
species_group_order = utils.species_group_order

# input
enrichment_file = "../data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/unique_term_count_within_databases.tsv"

# output
median_file               = "../data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/median_unique_term_count_within_databases/all_organisms.tsv"
mmu_median_file           = "../data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/median_unique_term_count_within_databases/10090.tsv"
hsa_median_file           = "../data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/median_unique_term_count_within_databases/9606.tsv"
nonhsa_nonmmu_median_file = "../data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/median_unique_term_count_within_databases/all_except_10090_9606.tsv"


enrichment_results = pd.read_table(enrichment_file)
enrichment_results.loc[:,"database"] = pd.Categorical(enrichment_results.database, 
                                                      categories = db_order,
                                                      ordered = True)
enrichment_results.loc[:,"species_group"] = pd.Categorical(enrichment_results.species_group,
                                                           categories = species_group_order,
                                                           ordered = True)

groups = enrichment_results.groupby(["species_group","database"])
median_results = groups.uniqueTermCount.apply(np.median).rename("median_term_count_per_user_input").reset_index()

median_results.to_csv(median_file, sep = "\t", index = False)


# median over human
hsa_results = enrichment_results.loc[enrichment_results.species_group.isin(["H.sapiens"]),:]
groups = hsa_results.groupby(["database"])
median_hsa = groups.uniqueTermCount.apply(np.median).rename("median_term_count_per_user_input").reset_index()
median_hsa.to_csv(hsa_median_file, sep = "\t", index = False)


# median over mouse
mmu_results = enrichment_results.loc[enrichment_results.species_group.isin(["M.musculus"]),:]
groups = mmu_results.groupby(["database"])
median_mmu = groups.uniqueTermCount.apply(np.median).rename("median_term_count_per_user_input").reset_index()
median_mmu.to_csv(mmu_median_file, sep = "\t", index = False)


# median over all organisms except mouse and human
nonhsa_nonmmu_results = enrichment_results.loc[~enrichment_results.species_group.isin(["H.sapiens", "M.musculus"]),:]
groups = nonhsa_nonmmu_results.groupby(["database"])
median_nonhsa_nonmmu_results = groups.uniqueTermCount.apply(np.median).rename("median_term_count_per_user_input").reset_index()
median_nonhsa_nonmmu_results.to_csv(nonhsa_nonmmu_median_file, sep = "\t", index = False)



