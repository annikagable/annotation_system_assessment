import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import unique_pathways as up
import utils
import sys

_, file_sigTermDf, file_nonredundant_term_count = sys.argv

# file_sigTermDf = "data/results/cameraPR/aggregation/sigTermDf.tsv"
# file_nonredundant_term_count = "data/results/cameraPR/redundancy_and_novelty/unique_term_count_within_databases.tsv"

dbColors = utils.dbColors
databases = dbColors.index


## import significantly enriched terms
sigTermDf = pd.read_table(file_sigTermDf)

## preprocessing

#sigTermDf = sigTermDf.loc[~sigTermDf.taxId.isin([9606, 10090]),:]

genes_series = sigTermDf.overlap_genes.apply(lambda string: set(string.split()))
sigTermDf = sigTermDf.assign( genes = genes_series).drop(columns = 'overlap_genes')
sigTermDf.rename(columns = {'overlap': 'nr_overlap_genes'}, inplace = True)


## remove redundant terms within own database
sigTermDf_only_unique_within_db = up.get_unique_within_database(sigTermDf = sigTermDf,
                                                                databases = databases)

## count number of nonredundant enriched terms within each database, per user input and species group
unique_terms_species_df = up.count_unique_terms(sigTermDf_unique_within_db = sigTermDf_only_unique_within_db,
                                                sigTermDf = sigTermDf,
                                                isUniqueAcross = False,
                                                databases = databases)

unique_terms_species_df.to_csv(file_nonredundant_term_count, sep = '\t', index = False)










