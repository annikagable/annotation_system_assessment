# Plot enriched terms per user input for each organism and database
# Terms which are subsets of other terms within the same database are omitted

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import utils
import sys

dbColors = utils.dbColors
sns.set_context('talk')

_,  file_nonredundant_term_count,  plot_file, plot_reduced_file = sys.argv 

# file_nonredundant_term_count = "data/results/cameraPR/redundancy_and_novelty/unique_term_count_within_databases.tsv"
# plot_file         = "figures/cameraPR/redundancy_and_novelty/all_species/nr_sig_terms_per_user_input.svg"
# plot_reduced_file = "figures/cameraPR/redundancy_and_novelty/reduced_species/nr_sig_terms_per_user_input.svg"


unique_term_counts = pd.read_table(file_nonredundant_term_count)

ylabel = "Non-redundant enrichment results – \n# of terms detected as significant per user input"
########################
## plot all organisms ##
########################

g = plt.figure(figsize = (25,10))
g = sns.boxplot(data = unique_term_counts,
                x = 'species_group',
                y = 'uniqueTermCount',
                hue = "database",
                palette = dbColors, 
                order = utils.species_group_order);
plt.yscale('symlog')
g.set_xlabel('')
g.set_ylabel(ylabel)
handles, labels = g.axes.get_legend_handles_labels()
g.legend_.remove()
g.legend(handles, labels, ncol=4, loc='upper right', frameon=True);
utils.savefig_multiformat(filename = plot_file)

#####################################
## plot a reduced set of organisms ##
#####################################

## Remove some species
unique_term_counts_reduced = unique_term_counts.loc[~unique_term_counts.species_group.isin(utils.drop_species),:].copy()

## Rename "other" category
unique_term_counts_reduced.loc[unique_term_counts_reduced.species_group == 'other', 'species_group'] = utils.rename_category['other']

## Plot
g = plt.figure(figsize = (20,10))
bp = sns.boxplot(data = unique_term_counts_reduced,
                x = 'species_group',
                y = 'uniqueTermCount',
                hue = "database",
                palette = dbColors, 
                order = utils.species_group_order_reduced,
                #linewidth = 1.2,
                fliersize = 3
                );
plt.yscale('symlog')
bp.set_xlabel('')
bp.set_ylabel(ylabel)
handles, labels = bp.axes.get_legend_handles_labels()
bp.legend_.remove()
bp.legend(handles, labels, ncol=4, loc='upper right', frameon=True);
utils.savefig_multiformat(filename = plot_reduced_file)
