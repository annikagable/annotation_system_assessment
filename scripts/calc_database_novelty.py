import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import unique_pathways as up
import utils
import sys

_, file_sigTermDf, file_sigTermDf_novelty, file_novel_term_count = sys.argv   

# file_sigTermDf = "data/results/cameraPR/aggregation/sigTermDf.tsv"
# out_dir = "data/results/cameraPR/redundancy_and_novelty" 
# file_sigTermDf_novelty = os.path.join(out_dir, 'sigTermDf_unique_terms_across_and_within_databases.tsv')
# file_novel_term_count  = os.path.join(out_dir, 'unique_term_count_across_and_within_databases.tsv'

dbColors = utils.dbColors
databases = dbColors.index

## import significantly enriched terms
sigTermDf = pd.read_table(file_sigTermDf)

## select the desired organism or data size
# sigTermDf = sigTermDf.loc[~sigTermDf.taxId.isin([9606, 10090]),:]

## convert the genes column into sets instead of string
genes_series = sigTermDf.overlap_genes.apply(lambda string: set(string.split()))
sigTermDf = sigTermDf.assign( genes = genes_series).drop(columns = 'overlap_genes')
sigTermDf.rename(columns = {'overlap': 'nr_overlap_genes'}, inplace = True)


## get enriched pathways unique across databases, for each user input separately
group = sigTermDf.groupby('dataId', sort = False, as_index = False)
sigTermDf_unique_terms = group.apply(up.determine_uniqueness_across_databases_per_term, 
                                     databases = databases).reset_index(drop = True)


## remove terms not unique across databases
sigTermDf_only_unique_terms = sigTermDf_unique_terms.loc[sigTermDf_unique_terms.isUniqueTerm, :].copy()


## remove redundant terms within own database
sigTermDf_unique_terms_within_database = up.get_unique_within_database(sigTermDf = sigTermDf_only_unique_terms, databases = databases)

## save to file
sigTermDf_unique_terms_within_database.to_csv(file_sigTermDf_novelty,
                                              sep = '\t',
                                              index = False)

## in case this dataframe is loaded from memory, it is important to make
## the database column categorical before proceeding:
# sigTermDf_unique_terms_within_database = pd.read_table(out_file)
# sigTermDf_unique_terms_within_database.loc[:,'database'] = pd.Categorical(sigTermDf_unique_terms_within_database.database, 
#                                                                           ordered = True,
#                                                                           categories = dbColors.index)


## count the number of unique terms per database and dataId
unique_terms_species_df = up.count_unique_terms(sigTermDf_unique_within_db = sigTermDf_unique_terms_within_database,
                                             sigTermDf = sigTermDf,
                                             isUniqueAcross = True,
                                             databases = databases)

unique_terms_species_df.to_csv(file_novel_term_count,
                               sep = '\t',
                               index = False)


## plot
## unfortunately, this includes so many zeros that the result looks very crappy
# up.plot_term_counts(data = unique_terms_species_df,
#                      palette = dbColors,
#                      ylabel = 'added novelty\n(number of terms unique across databases)',
#                      plot_file_tag = "unique_terms_novelty",
#                      out_dir = out_dir)
# 
# 
# ## alternative, non-logged plots that are able to (sort of) display the zeros:
# 
# data = unique_terms_species_df
# palette = dbColors
# ylabel = 'added novelty\n(number of terms unique across databases)'
# plot_file_tag = "unique_terms_novelty"
# 
# 
# # violin
# plt.figure(figsize = (6,4))
# sns.violinplot(data = data, 
#                x = 'database', 
#                y = 'uniqueTermCount', 
#                palette = palette,
#                inner = "box", #{"box", "quartile", "point", "stick", None}
#                scale = "count",#{"area", "count", "width"}
#                cut = 0);
# plt.ylim(-0.5,5)
# plt.xticks(rotation=45, horizontalalignment='right')
# plt.xlabel('');
# plt.ylabel(ylabel);
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, f'{plot_file_tag}.violin.png'))
# 
# 
# 
# # box
# # PubMed: up to 18000 unique terms added per dataset
# plt.figure(figsize = (6,4))
# sns.boxplot(data = data, 
#                x = 'database', 
#                y = 'uniqueTermCount', 
#                palette = palette,
#                fliersize=1);
# plt.ylim(-0.5,40)
# plt.xticks(rotation=45, horizontalalignment='right')
# plt.xlabel('');
# plt.ylabel(ylabel);
# plt.text(x = 10, y = 42, s = '*', fontsize='large', horizontalalignment='center')
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, f'{plot_file_tag}.box.png'))
# 
# 
# 
# # boxen
# plt.figure(figsize = (6,4))
# sns.boxenplot(data = data, 
#                x = 'database', 
#                y = 'uniqueTermCount', 
#                palette = palette);
# plt.ylim(-0.5,40)
# plt.xticks(rotation=45, horizontalalignment='right')
# plt.xlabel('');
# plt.ylabel(ylabel);
# plt.text(x = 10, y = 42, s = '*', fontsize='large', horizontalalignment='center')
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, f'{plot_file_tag}.boxen.png'))
# 
# 
# 
# 
# 
# ## count the number of dataIds which have any unique enrichments
# 
# #unique_term_count_sum = unique_terms_species_df.groupby('dataId').uniqueTermCount.sum()
# #len(unique_term_count_sum)
# #len(unique_term_count_sum[unique_term_count_sum != 0])
# 
# 
# ## plot "any added novelty"
# ## aka does the database provide at least one unique term for the given input? 
# u = unique_terms_species_df.assign(isUnique = unique_terms_species_df.uniqueTermCount > 0)
# sns.countplot(data = u, hue = 'isUnique', y = 'database', palette = ['tab:red', 'dodgerblue']);
# plt.ylabel("");
# plt.xlabel("# of user input datasets");
# plt.legend(title= "at least one\nunique term")
# plt.tight_layout()
# plot_file_tag = "unique_terms_count_any_novelty"
# plt.savefig(os.path.join(out_dir, f'{plot_file_tag}.png'))
# 

