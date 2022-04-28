import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import unique_pathways as up
import utils
import sys

_, file_sigTermDf, file_sigTermDf_novelty, file_novel_term_count = sys.argv   

# file_sigTermDf = "data/results/cameraPR/overlap_3-200/aggregation/sigTermDf_alpha0.05.tsv"
# out_dir = "data/results/cameraPR/overlap_3-200/redundancy_and_novelty" 
# file_sigTermDf_novelty = os.path.join(out_dir, 'sigTermDf_unique_terms_across_and_within_databases.tsv')
# file_novel_term_count  = os.path.join(out_dir, 'unique_term_count_across_and_within_databases.tsv')

dbColors = utils.dbColors
databases = dbColors.index

unique_within_db_columns = ['term_shorthand', 
                             'dataId',
                             'database', 
                             'taxId', 
                             'species',
                             'species_group', 
                             'nr_overlap_genes', 
                             'nr_unique_genes', 
                             'isUniqueTerm',
                             'genes', 
                             'unique_genes']


## import significantly enriched terms
sigTermDf = pd.read_table(file_sigTermDf)


if len(sigTermDf) == 0:
    # Prevent any errors due to empty input dataframe
    sigTermDf_unique_terms_within_database = pd.DataFrame(columns = unique_within_db_columns)
    unique_terms_species_df                = pd.DataFrame(columns = up.unique_per_db_and_species_column_order)


else:

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


    
    ## remove redundant terms within own database (discard enriched term subsets)
    sigTermDf_unique_terms_within_database = up.get_unique_within_database(sigTermDf = sigTermDf_only_unique_terms,
                                                                           databases = databases)
    assert (sigTermDf_unique_terms_within_database.columns == unique_within_db_columns).all()
    
    

    ## count the number of unique terms per database and dataId
    unique_terms_species_df = up.count_unique_terms(sigTermDf_unique_within_db = sigTermDf_unique_terms_within_database,
                                                     sigTermDf = sigTermDf,
                                                     isUniqueAcross = True,
                                                     databases = databases,
                                                   column_order = up.unique_per_db_and_species_column_order)

    

## save to file
sigTermDf_unique_terms_within_database.to_csv(file_sigTermDf_novelty,
                                              sep = '\t',
                                              index = False)

unique_terms_species_df.to_csv(file_novel_term_count,
                               sep = '\t',
                               index = False)



