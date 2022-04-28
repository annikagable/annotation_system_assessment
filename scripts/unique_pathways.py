import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

unique_per_db_and_species_column_order = ['dataId', 
                                          'database', 
                                          'uniqueTermCount', 
                                          'taxId', 
                                          'species',
                                          'species_group']
    
def discard_term_subsets(df):
    unique_terms = []
    isUnique_list = []
    for row in df.itertuples():
        isExisting = False
        for t in unique_terms:
            if row.genes.issubset(t):
                isExisting = True
                break
        if not isExisting:
            unique_terms.append(row.genes)
        isUnique_list.append(not isExisting)

    filtered_df = df.loc[isUnique_list,:].copy()
    return filtered_df
    
def get_unique_within_database(sigTermDf, databases):
    '''
    sigTermDf: A dataframe of significantly enriched terms with at least the columns 
                'term_shorthand', 'dataId', 'taxId', 'species', 'species_group', 
                'database', 'nr_overlap_genes'.
    '''
    
    
    sigTermDf_sorted = sigTermDf.sort_values(by = ['dataId',
                                                   'taxId', 
                                                   'species', 
                                                   'species_group', 
                                                   'database',
                                                   'nr_overlap_genes'],
                                                    ascending = False)
    
    sigTermDf_unique = sigTermDf_sorted.groupby(['dataId',
                                                 'taxId', 
                                                 'species', 
                                                 'species_group', 
                                                 'database'],
                                  as_index = False).apply(discard_term_subsets).reset_index(drop=True)
    
    sigTermDf_unique.loc[:, 'database'] = pd.Categorical(sigTermDf_unique.database,
                                                           ordered = True,
                                                           categories = databases)
    
    
    return sigTermDf_unique




def determine_uniqueness_across_databases_per_term(sigTermDf, databases, min_unique_ratio = 0.5 ):
    
    '''
    get enriched pathways unique across databases, for each user input separately
    '''
    
    
    #import pdb; pdb.set_trace()
    records = []
    for row in sigTermDf.itertuples():

        # Get the union of all genes enriched for this dataId in all the terms from the other databases
        other_genes_list = sigTermDf.loc[sigTermDf.database != row.database, 'genes'].tolist()
        other_genes = set().union(*other_genes_list)

        # Get the genes which are unique to this term
        unique_genes = row.genes - other_genes
        nr_unique_genes = len(unique_genes)

        # Is the ratio of unique genes to total genes in the term more than half?
        if nr_unique_genes/row.nr_overlap_genes >= min_unique_ratio:
            isUniqueTerm = True
        else:
            isUniqueTerm = False

        # Create a dataframe
        records.append((row.term_shorthand,
                        row.dataId,
                        row.database,
                        row.taxId,
                        row.species,
                        row.species_group,
                        row.nr_overlap_genes,
                        nr_unique_genes,
                        isUniqueTerm,
                        row.genes,
                        unique_genes
                        ))

    df_columns = ['term_shorthand', 'dataId', 'database', 'taxId', 'species', 'species_group', 
                  'nr_overlap_genes', 'nr_unique_genes', 'isUniqueTerm', 'genes', 'unique_genes']
    sigTermDf_unique_terms = pd.DataFrame.from_records(records, columns = df_columns)
    sigTermDf_unique_terms.loc[:, 'database'] = pd.Categorical(sigTermDf_unique_terms.database,
                                                               ordered = True,
                                                               categories = databases)
    return sigTermDf_unique_terms





    

def count_unique_terms(sigTermDf_unique_within_db, sigTermDf, isUniqueAcross, databases, column_order):
        
    if isUniqueAcross:
        ## count the number of unique terms for each dataId and each database
        unique_terms_df = sigTermDf_unique_within_db.groupby(['dataId', 'database'], 
                                                              as_index = False,
                                                              observed = False).isUniqueTerm.sum()
        unique_terms_df.rename(columns = {'isUniqueTerm': 'uniqueTermCount'}, inplace = True)
        
    else:
        
        important_columns = ['dataId','taxId', 'species_group', 'database']
        unique_terms_df = sigTermDf_unique_within_db.groupby(important_columns,
                                                 as_index = False,
                                                 observed = False).q_value.count()
        unique_terms_df.rename(columns = {'q_value': 'uniqueTermCount'}, inplace = True)
    
    
    # Set databases with no enriched terms at all to 0 instead of NaN
    unique_terms_df.loc[unique_terms_df.uniqueTermCount.isna(), 'uniqueTermCount'] = 0



    ## add the totally non-unique dataIds back in which dropped out
    all_dataIds = set(sigTermDf.dataId.unique())
    unique_dataIds = set(unique_terms_df.dataId.unique())
    dropped_dataIds = all_dataIds - unique_dataIds
    records = []
    for d in dropped_dataIds:
        for db in databases:
            records.append((d,db,0))
    missing_dataId_df = pd.DataFrame.from_records(records, 
                                                  columns = ['dataId', 
                                                             'database',
                                                             'uniqueTermCount'])
    unique_terms_df = unique_terms_df.append(missing_dataId_df)


    ## add the species and taxId back in 
    dataId_to_species = sigTermDf.loc[:, ['dataId', 'taxId', 'species', 'species_group']].drop_duplicates()
    unique_terms_species_df = unique_terms_df.merge(dataId_to_species)
    
    unique_terms_species_df.loc[:,'database'] = pd.Categorical(unique_terms_species_df.database, 
                                                              ordered = True,
                                                              categories = databases)
    
    ## sort the columns
    unique_terms_species_df = unique_terms_species_df.loc[:, column_order]
    
    
    return unique_terms_species_df






def plot_term_counts(data,
                     palette,
                     ylabel,
                     plot_file_tag,
                     out_dir):
    '''
    data:            unique_terms_species_df, the dataframe with the uniqueTermCount column
    ylabel:          use 'added novelty\n(number of terms unique across databases)' or
                         '# of unique significant terms per user input'
    plot_file_tag:   'unique_terms_novelty' or 'unique_terms_within'
    out_dir:         the output dir for the plots
    '''

    plotting_funs = {'box':   sns.boxplot,
                     'boxen': sns.boxenplot}
    
    for fun in plotting_funs:
        plt.figure(figsize = (6,4))
        plotting_funs[fun](data = data, x = 'database', y = 'uniqueTermCount', palette = palette);
        plt.yscale('symlog');
        plt.xticks(rotation=45, horizontalalignment='right')
        plt.xlabel('');
        plt.ylabel(ylabel);
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f'{plot_file_tag}.{fun}.log.png'))



    
    
