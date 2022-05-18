import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import utils


_, sigTermDf_file, species_members_file, taxId, cum_gene_unions_file, summary_cum_gene_unions_file = sys.argv

# taxId = 9606

# ## Input files
# sigTermDf_file = "data/results/cameraPR/overlap_3-200/aggregation/sigTermDf_alpha0.05.tsv" 
# species_members_file = "data/raw/global_enrichment_annotations/9606.terms_members.tsv"

## Output files
# # computed on all possible term sizes => for plotting summary curves
# summary_cum_gene_unions_file = "data/results/cameraPR/overlap_3-200/unique_enriched_genes/9606/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv"
# # computed only on existing term sizes => for plotting individual users curves
# cum_gene_unions_file= "data/results/cameraPR/overlap_3-200/unique_enriched_genes/9606/cumulative_gene_unions_by_term_size.tsv"


# sigTermDf_file = "data/results/cameraPR/overlap_0-Inf/aggregation/sigTermDf_alpha0.05.tsv" 
# species_members_file = "data/raw/global_enrichment_annotations/9606.terms_members.tsv"
# summary_cum_gene_unions_file = "data/results/cameraPR/overlap_0-Inf/unique_enriched_genes/9606/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv"
# cum_gene_unions_file= "data/results/cameraPR/overlap_0-Inf/unique_enriched_genes/9606/cumulative_gene_unions_by_term_size.tsv"


def uniquify(df_term_sizes):

    # get the unique genes by term size, dataId, database
    grouping_columns = ['term_size', 'species', 'species_group', 'database', 'dataId']
    group = df_term_sizes.groupby(grouping_columns, as_index = False)
    gene_unions_by_term_size = group.genes.apply(lambda gs: set.union(*gs)).reset_index(name = "gene_union")


    # get the number of unique genes
    gene_union_size = gene_unions_by_term_size.gene_union.apply(len)
    gene_unions_by_term_size = gene_unions_by_term_size.assign(gene_union_size = gene_union_size)
    
    return gene_unions_by_term_size


def get_cumulative_union(df):

    # sort by the term size
    df = df.sort_values('term_size')
    
    # initialize: compute gene union for first term size
    unions = [np.nan]*len(df)
    t0 = df.term_size.values[0]
    gene_sets = df.loc[(df.term_size == t0), 'gene_union'].to_list()
    unions[0] = set.union(*gene_sets)
    
    for i, t in list(enumerate(df.term_size))[1:]:
        # compute gene union based on previous unions
        gene_sets = df.loc[(df.term_size == t), 'gene_union'].to_list()
        gene_set_list = gene_sets + [unions[i-1]] 
        unions[i] = set.union(*gene_set_list)
        
    union_sizes = [len(u) for u in unions]
    
    return df.assign(cumulative_gene_union = unions,
                     cumulative_gene_union_size = union_sizes)
    
    
def uniquify_thresholded(df_term_sizes):
    
    # get the unique cumulative genes by term_size, dataId and database
    grouping_columns = ['species', 'species_group', 'database', 'dataId']
    group = df_term_sizes.groupby(grouping_columns, group_keys = False, as_index = False)
    gene_unions_by_term_size = group.apply(get_cumulative_union) # alternative, same result
    
    return gene_unions_by_term_size


def get_percentiles(df, size_column, grouping_columns, percentiles):
    
    sizes = df.loc[:,size_column]
    l, u = percentiles
    
    lower = np.percentile(sizes, l)
    upper = np.percentile(sizes, u)
    med = np.median(sizes)
    av = np.mean(sizes)
    
    # make sure that there are no decimal dots in the column names
    l = str(l).replace('.', '_')
    u = str(u).replace('.', '_')
        
    df = pd.DataFrame({'term_size': df.term_size.iat[0],
                       'species': df.species.iat[0],
                       'species_group': df.species_group.iat[0],
                       'database': df.database.iat[0],
                       f'perc_{l}_gene_union_size': [lower], 
                       f'perc_{u}_gene_union_size': [upper],
                       'median_gene_union_size': [med],
                       'mean_gene_union_size': [av]})

    return df




dbColors = utils.dbColors
percentile_band = (47.5, 52.5)
taxId = int(taxId)

# define which columns the output dataframes should have
gene_unions_by_term_size_columns = ['term_size', 'species', 'species_group', 
                                    'database', 'dataId',
                                    'gene_union', 'gene_union_size']


l, u = [str(p).replace('.', '_') for p in percentile_band]
summary_cum_gene_union_df_columns =['term_size', 'species', 'species_group', 
                                    'database',
                                   f'perc_{l}_gene_union_size', f'perc_{u}_gene_union_size',
                                   'median_gene_union_size', 'mean_gene_union_size']

cum_gene_unions_by_term_size_columns = gene_unions_by_term_size_columns + ['cumulative_gene_union',
                                                                           'cumulative_gene_union_size']

# import enriched terms
sigTermDf = pd.read_table(sigTermDf_file)
sigTermDf.loc[:,'database'] = pd.Categorical(sigTermDf.database, 
                                             ordered = True,
                                             categories = dbColors.index)

# get terms of one species only
df = sigTermDf.loc[sigTermDf.taxId == taxId, :]


if len(df) == 0:
    gene_unions_by_term_size = pd.DataFrame(columns = gene_unions_by_term_size_columns)

else:
    
    species       = df.species.iloc[0]       # e.g. 'Homo sapiens'
    species_group = df.species_group.iloc[0] # e.g. 'H.sapiens'

    ## convert the genes column into sets instead of string
    genes_series = df.overlap_genes.apply(lambda string: set(string.split()))
    df = df.assign( genes = genes_series).drop(columns = 'overlap_genes')
    df.rename(columns = {'overlap': 'nr_overlap_genes'}, inplace = True)


    ## import raw term sizes
    species_members = pd.read_table(species_members_file, usecols = [0,2], names = ["term_shorthand", "term_size"])

    # join enriched terms with raw term sizes
    df_term_sizes = df.merge(species_members)

    # get enriched gene unions by database and term size
    gene_unions_by_term_size = uniquify(df_term_sizes)
    # => this table is needed both for calculations for all possible term sizes (summary plots)
    # and for calculations on only existing term sizes (individual users' plots)

    assert np.all(gene_unions_by_term_size.columns == gene_unions_by_term_size_columns)

# # save to file
# gene_unions_by_term_size.to_csv(gene_unions_by_term_size_file, 
#                               sep = "\t", 
#                               index = False)



###########################################################################################################
######## Get gene union size only for the term sizes that exist in each user input and database ###########
###########################################################################################################

# In order to plot individual user's curves, we only want to include the term sizes that actually appear in the
# enrichment

if len(df) == 0:
    # create empty dataframe
    cum_gene_unions_by_term_size = pd.DataFrame(columns = cum_gene_unions_by_term_size_columns)
else:
    cum_gene_unions_by_term_size = uniquify_thresholded(gene_unions_by_term_size)
    assert np.all(cum_gene_unions_by_term_size.columns == cum_gene_unions_by_term_size_columns)

# save to file
cum_gene_unions_by_term_size.to_csv(cum_gene_unions_file, 
                                      sep = "\t", 
                                      index = False)



####################################################################################################
######## Get the summary of all user inputs: calc. union size at all possible term sizes ###########
####################################################################################################

# Not all dataIds have enrichments for any given term size and database.
# Therefore I need to add missing combinations of term size and database

if len(df) == 0:
    # create empty dataframe
    summary_cum_gene_union_df = pd.DataFrame(columns = summary_cum_gene_union_df_columns)
else:
    
    # Get all possible combinations of term size / database / dataId
    all_databases = dbColors.index.to_list()
    all_dataIds = df.dataId.unique()

    records = []
    for db in all_databases:
        all_term_sizes = gene_unions_by_term_size.loc[(gene_unions_by_term_size.database == db),'term_size'].unique()
        for term_size in all_term_sizes:
            for dataId in all_dataIds:
                    records.append((term_size, db, dataId))

    combos = pd.DataFrame.from_records(records, columns = ['term_size', 'database', 'dataId'])

    # Join combos with existing data
    gene_unions_by_term_size_added_missing = gene_unions_by_term_size.merge(combos, how = 'outer')

    # Replace NaNs generated by outer join with the proper values
    na_condition = gene_unions_by_term_size_added_missing.gene_union.isna()
    constant_column_names = ['species', 'species_group', 'gene_union', 'gene_union_size']
    constant_column_values = [species, species_group, set(), 0]
    gene_unions_by_term_size_added_missing.loc[na_condition, constant_column_names] = constant_column_values

    # Get the cumulative union sizes on the data where missing combos were addded
    cum_gene_unions_by_all_term_sizes = uniquify_thresholded(gene_unions_by_term_size_added_missing)


    # Get the average by term size and database for easier/faster plotting
    grouping_columns = ['term_size', 'species', 'species_group', 'database']
    cgubts_group = cum_gene_unions_by_all_term_sizes.groupby(grouping_columns, as_index = False)
    summary_cum_gene_union_df = cgubts_group.apply(get_percentiles,
                                                size_column = 'cumulative_gene_union_size',
                                                grouping_columns = grouping_columns,
                                                percentiles = percentile_band,
                                                )
# Save to file
summary_cum_gene_union_df.to_csv(summary_cum_gene_unions_file, 
                              sep = "\t", 
                              index = False)


    

