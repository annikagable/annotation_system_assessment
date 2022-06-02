import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import utils
import sys

def plot_novelty(any_unique_terms_data, input_counts, ylabel):
    
    h = sns.catplot(data = any_unique_terms_data,
                kind = "count",
                x = "database",
                col = "species_group",
                palette = dbColors,
                sharey=False,
                height = 5,
                aspect = 0.4);
    
    h.set_titles('')#col_template = '{col_name}')
    
    
    for i, ax in enumerate(h.axes[0]):
        ax.axhline(input_counts[i], ls='--', label="total # of user inputs")
        ax.set_xlabel(input_counts.index[i], labelpad = 18)
        ax.set_xticklabels('')
        ax.set_xticks([])
        ax.yaxis.get_major_locator().set_params(integer=True)
        
    h.axes[0,0].set_ylabel(ylabel)
    


_, dataId_isSig_file, unique_sigTermDf_file, novelty_plot_file, novelty_plot_reduced_file, novel_count_file = sys.argv

# dataId_isSig_file     = "data/results/cameraPR/overlap_3-Inf/aggregation/dataId_isSig_alpha0.05.tsv"
# unique_sigTermDf_file = "data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/sigTermDf_unique_terms_across_and_within_databases.tsv"
# novelty_plot_file         = "figures/cameraPR/overlap_3-Inf/redundancy_and_novelty/all_species/at_least_one_novel_significant.svg"
# novelty_plot_reduced_file = "figures/cameraPR/overlap_3-Inf/redundancy_and_novelty/reduced_species/at_least_one_novel_significant.svg"
# novel_count_file = "data/results/cameraPR/overlap_3-Inf/redundancy_and_novelty/novel_enriched_count.tsv"

sns.set_context('talk')

dbColors = utils.dbColors
species_group_order = utils.species_group_order
drop_species = utils.drop_species
rename_category = utils.rename_category

## Count how many datasets have at least one novel term, per species_group and database ##

# Load the table with novel significant terms
sigTermDf_unique = pd.read_table(unique_sigTermDf_file)
sigTermDf_unique.loc[:,'database'] = pd.Categorical(sigTermDf_unique.database, 
                                                    categories = dbColors.keys(), 
                                                    ordered = True)
sigTermDf_unique.loc[:,'species_group'] = pd.Categorical(sigTermDf_unique.species_group, 
                                                    categories = species_group_order, 
                                                    ordered = True)

if len(sigTermDf_unique) == 0:
    any_unique_terms_by_database = pd.DataFrame(columns = ['species_group', 'database', 'dataId', 'term_shorthand'])
    any_unique_terms_by_database.loc[:,'database'] = pd.Categorical(any_unique_terms_by_database.database, 
                                                        categories = dbColors.keys(), 
                                                        ordered = True)
    any_unique_terms_by_database.loc[:,'species_group'] = pd.Categorical(any_unique_terms_by_database.species_group, 
                                                        categories = species_group_order, 
                                                        ordered = True)
    
else:
    # Count datasets with any unique terms
    group = sigTermDf_unique.groupby(["species_group","database", "dataId"])
    count_any_unique = group.term_shorthand.count()
    any_unique_terms_by_database = (count_any_unique > 0)
    any_unique_terms_by_database = any_unique_terms_by_database[any_unique_terms_by_database].reset_index()


novel_enriched_input_count = any_unique_terms_by_database.groupby(["species_group", "database"]).apply(len).rename("novel_enriched_count").reset_index()
novel_enriched_input_count.to_csv(novel_count_file, sep = "\t", index = False)


## Get total number of inputs in each species group ##
dataId_isSig = pd.read_table(dataId_isSig_file)
dataId_isSig.loc[:,'species_group'] = pd.Categorical(dataId_isSig.species_group,
                                                    categories = species_group_order, 
                                                    ordered = True)
input_counts = dataId_isSig.groupby("species_group").dataId.unique().apply(len)
input_counts.name = "all_input_counts"



## Remove two species and rename "other" species category ##
any_unique_terms_by_database_reduced = utils.remove_species_and_rename_other(df = any_unique_terms_by_database,
                                                                             drop_species = drop_species,
                                                                             rename_category = rename_category)

input_counts_reduced = input_counts[~input_counts.index.isin(drop_species)]
input_counts_reduced.rename(rename_category, inplace = True)



## Plot the datasets with at least one novel term ##
ylabel = """Added novelty – # of user inputs
 with terms not detected elsewhere"""

## all species
plot_novelty(any_unique_terms_data = any_unique_terms_by_database, 
                 input_counts = input_counts, 
                 ylabel = ylabel)
utils.savefig_multiformat(filename = novelty_plot_file)

## reduced species
plot_novelty(any_unique_terms_data = any_unique_terms_by_database_reduced, 
                 input_counts = input_counts_reduced, 
                 ylabel = ylabel)
utils.savefig_multiformat(filename = novelty_plot_reduced_file)
