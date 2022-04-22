# This script reads the enrichment results from both limma's cameraPR and Global Enrichment. 
# In cameraPR, each database was a separate enrichment job, whereas in Global Enrichment, all databases are in
# one combined output file. 
# All significant results from all species, all dataIds, all databases, are combined into one big dataframe.
# A separate dataframe is generated to keep all species and dataIds where there was no significant result at all,
# so that we can compare which percentage of the user inputs got a significant enrichment.

import os
import sys
import pandas as pd
import numpy as np
import utils
    
    
    
def read_one_user_output_cameraPR(f, alpha):
    
    df = pd.read_table(f)
    
    if df.columns[0].startswith('data'):
        # Catch a common error when manually updating the aggregated data.
        print(f"A wrong file is listed in the aggregated file: {f}")
        sys.exit()
            
#     if "q_value" not in df.columns:
#         ## This can happen if there is only one term.
#         df = df.assign(q_value = df.p_value)

        
    df = df.loc[df.q_value < alpha,:]
        
    ## get unique identifier & species ID of the user input
    basename = os.path.basename(f)
    dataId, taxId, database, _ = basename.split('.')

    df = df.assign(dataId = dataId,
                   taxId = int(taxId),
                   database = database)
    
        
    # also build a list of all dataIds with their taxIds separately, so that 
    # we can check which databases are missing which dataIds/species
    atLeastOneSigTerm = False if len(df) == 0 else True
    dataDf = pd.DataFrame({"dataId": dataId,
                           "taxId": int(taxId),
                           "database": database,
                           "atLeastOneSigTerm": atLeastOneSigTerm}, index = [0])
    return df, dataDf



def read_one_user_output_GE(f, alpha):

    column_names = ['method',  'direction', 'enrichment_score',
                    'term',    'etype',     'overlap',
                    'p_value', 'q_value',   'shorthands']
    
    # prevent random dtype assignments 
    column_dtypes = [np.dtype('O'),
                     np.dtype('int64'),
                     np.dtype('float64'),
                     np.dtype('int64'),
                     np.dtype('int64'),
                     np.dtype('int64'),
                     np.dtype('float64'),
                     np.dtype('float64'),
                     np.dtype('O')]

    col_dtype_dict = dict(zip(column_names, column_dtypes))
    
    df = pd.read_table(f, header = None, usecols=range(9),
                       names = column_names, dtype = col_dtype_dict)

    df = df.loc[df.q_value < alpha,:]
        
    ## get unique identifier & species ID of the user input
    basename = os.path.basename(f)
    dataId = basename.split('.')[0]
    taxId = int(basename.split('.')[1])
    
    df = df.assign(dataId = dataId, taxId = taxId)
    
    ## add database
    df = utils.add_database_column(df, col= "etype", translate_to = 'database',
                                   etype_to_database = utils.etype_to_database,
                                   order = utils.dbColors.index)
    
    
    # also build a list of all dataIds with their taxIds separately, so that 
    # we can check which databases are missing which dataIds/species
        
    termCount = df.groupby("database").term.apply(len)
    atLeastOneSigTerm = termCount > 0
    
    dataDf = pd.DataFrame({"dataId": dataId,
                           "taxId": int(taxId),
                           "database": db_order,
                           "atLeastOneSigTerm": atLeastOneSigTerm})

    return df, dataDf




def read_all_user_outputs(output_files_file,
                          enrichment_method, 
                          db_order,
                          alpha):
    
    
    with open(output_files_file, 'r') as files_f:
        output_files = [l.strip() for l in files_f.readlines()]

    
    if enrichment_method == "cameraPR":
        read_one_user_output = read_one_user_output_cameraPR
    elif enrichment_method == "GE":
        read_one_user_output = read_one_user_output_GE
    else:
        sys.exit("enrichment_method needs to be 'cameraPR' or 'GE'.")
        
    dfList, dataDfList = zip(*[read_one_user_output(f, alpha) for f in output_files])

    print("files imported.")

    termDf = pd.concat(dfList, ignore_index=True)
    
    dataIdDf = pd.concat(dataDfList, ignore_index=True)
    
    ## Keep only databases which are actually present
    db_order = np.array(db_order)
    index = np.isin(db_order, dataIdDf.database.unique())
    cat_order = db_order[index]
       
    dataIdDf.loc[:,'database'] = pd.Categorical(dataIdDf.database, categories = db_order, ordered = True)
    termDf.loc[:,'database'] = pd.Categorical(termDf.database, categories = db_order, ordered = True)
    
    return termDf, dataIdDf




def add_species_names(df, species_taxIds_file):
    """
    Add a "species" column (containing the full species name) to a dataframe that has a "taxId" column.
    
    termDf: the dataframe with the taxId column
    species_taxIds_file: a file containing a table with "taxId" and "species" columns
    
    Returns the merge of the two tables.
    """
    
    assert("taxId" in df.columns)
    
    column_names = ["taxId", "species"]
    column_dtypes = [int, str]
    col_dtype_dict = dict(zip(column_names, column_dtypes))

    speciesTaxIds = pd.read_table(species_taxIds_file, usecols= [0,2],
                                  names = column_names, header=0, dtype = col_dtype_dict)

    termNamedDf = df.merge(speciesTaxIds, how='left')
    
    return termNamedDf


def add_species_groups(df, n_grouped_species = 8):
    '''
    Use the n species with the most user inputs, and group all other
    species into an 'other' category. Add this species group column to
    the existing dataframe.
    '''
    
    assert(all([col in df.columns for col in ["dataId", "species", "taxId"]]))
    
    ### Part 1: Determine the most common species
    
    ## Get each dataId and the associated species
    dataAndSpecies = df.loc[:,["dataId", "species", "taxId"]].drop_duplicates()
    
    ## Count how many dataIds per species we have
    species_counts = dataAndSpecies.groupby(['species', 'taxId']).dataId.count().sort_values(ascending = False)

    
    ## Select the most frequent n species to be species groups
    species_group_long_names = species_counts.head(n_grouped_species).reset_index(level = 0).drop(columns = 'dataId')

    assert(type(species_group_long_names) == pd.DataFrame)
    
    
    ### Part 2: Create the groups from the most common species
    ### This part could also be applied to any dataframe as long
    ### as we know the species to choose.
    
    ## Create shortened names using the initial of the genus plus the species name
    def create_short_name(species_name):
        name_list = species_name.split()
        initial = name_list[0][0]
        species = name_list[1]
        return f"{initial}.{species}"

    species_group_short_names = species_group_long_names.species.apply(create_short_name)
    species_group_short_names.name = "species_group"
    species_group_dict  = species_group_short_names.to_dict()

    ## If the species has no group, it's group is other
    def species_to_group(tax):
        return species_group_dict.get(tax, "other")

    ## Add species group column
    df = df.assign(species_group = df.taxId.apply(species_to_group))
    
    ## Convert the species to categorical
    df.loc[:,'species_group'] = pd.Categorical(df.species_group, 
                                                    categories = species_group_short_names.to_list()+["other"], 
                                                    ordered = True)
    
    return df



# enrichment_files_file = "data/results/cameraPR/aggregation/aggregated.txt"
# species_taxIds_file = "data/raw/species.v11.0.txt"
# enrichment_method = "cameraPR"
# alpha = 0.05
# n_grouped_species = 8
# out_dir = "data/results/cameraPR/aggregation"

# enrichment_files_file = "data/results/GE/deduplicated_output_files.tsv"
# species_taxIds_file = "data/raw/species.v11.0.txt"
# enrichment_method = "GE"
# alpha = 0.05
# n_grouped_species = 8
# out_dir = "data/results/GE/aggregation"
                          

# # when using the aggregated data with downsampled PubMed for term size v effect size plot
# enrichment_files_file = "data/aggregation_downsampled_PubMed/aggregated.txt"
# species_taxIds_file = "../../reproducible_enrichment/data/raw/species.v11.0.txt"
# enrichment_method = "cameraPR"
# alpha = 1
# n_grouped_species = 8
# out_dir = "data/aggregation_downsampled_PubMed"


_, enrichment_files_file, species_taxIds_file, alpha, n_grouped_species, enrichment_method, out_dir = sys.argv
alpha = float(alpha)
n_grouped_species = int(n_grouped_species)

### Reading & preprocessing data

termDf, dataIdDf = read_all_user_outputs(output_files_file = enrichment_files_file,
                                         enrichment_method = enrichment_method,
                                         db_order = utils.dbColors.index,
                                         alpha = alpha)

termNamedDf   = add_species_names(df = termDf,   species_taxIds_file = species_taxIds_file)
dataIdNamedDf = add_species_names(df = dataIdDf, species_taxIds_file = species_taxIds_file)

termNamedGroupDf   = add_species_groups(df = termNamedDf,   n_grouped_species = n_grouped_species)
dataIdNamedGroupDf = add_species_groups(df = dataIdNamedDf, n_grouped_species = n_grouped_species)

## Put overlap gene columns to the end of the dataframe for easier display as tsv
overlap_gene_column = termNamedGroupDf.pop('overlap_genes')
termNamedGroupDf = pd.concat([termNamedGroupDf, overlap_gene_column], axis = 1)

sigTermDf_file    = os.path.join(out_dir, "sigTermDf_alpha"+str(alpha)+".tsv")
dataId_isSig_file = os.path.join(out_dir, "dataId_isSig_alpha"+str(alpha)+".tsv")

termNamedGroupDf.to_csv(sigTermDf_file, sep = '\t', index = False)
dataIdNamedGroupDf.to_csv(dataId_isSig_file, sep = '\t', index = False)
