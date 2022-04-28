import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial']

dbColors = pd.Series({"Reactome": "#FFA52B",
                      "KEGG": "#C79D36",
                      "UniProt": "#F0F01F",
                      "GO_BP": "#2ECC71",
                      "GO_CC": "#229954",
                      "GO_MF": "#44DEBD",
                      "SMART": "#AD86FF", 
                      "Pfam": "#8E44AD",
                      "InterPro": "#AF7AC5",
                      "STRINGclusters": "#E74C3C",
                      "PubMed": "#95A5A6"}, name = 'color')

dbDashes = {"Reactome": (8,1,2, 1),
              "KEGG": (5, 1),
              "UniProt": (1,0),
              "GO_BP": (2, 1),
              "GO_CC": (6, 1),
              "GO_MF": (3,1),
              "SMART": (4, 0.5, 1, 0.5), 
              "Pfam": (10, 1),
              "InterPro": (1,1),
              "STRINGclusters": (3, 1, 1, 1),
              "PubMed": (3, 2)}

dbNames = """Reactome pathways
KEGG pathways
UniProt keywords
GO Biological Process
GO Cellular Component
GO Molecular Function
SMART domains
Pfam families
InterPro entries
STRING clusters
tagged PubMed publications""".split("\n")

database_to_etype = { 
    'GO_BP':-21,
    'GO_CC':-22,
    'GO_MF':-23,
    'BTO':-25,
    'DOID':-26,
    'UniProt':-51,
    'KEGG':-52,
    'SMART':-53,
    'InterPro':-54,
    'Pfam':-55,
    'PubMed':-56,
    'Reactome':-57,
    'STRINGclusters':-78
}


etype_to_database = { 
    -78: 'STRINGclusters',
    -57: 'Reactome',
    -56: 'PubMed',
    -55: 'Pfam',
    -54: 'InterPro',
    -53: 'SMART',
    -52: 'KEGG',
    -51: 'UniProt',
    -23: 'GO_MF',
    -22: 'GO_CC',
    -21: 'GO_BP',
    ### 
    78: 'STRINGclusters',
    57: 'Reactome',
    56: 'PubMed',
    55: 'Pfam',
    54: 'InterPro',
    53: 'SMART',
    52: 'KEGG',
    51: 'UniProt',
    23: 'GO_MF',
    22: 'GO_CC',
    21: 'GO_BP'
}

species_group_order = np.array(['H.sapiens', 'M.musculus', 'R.norvegicus', 'D.melanogaster',
                       'A.thaliana','D.rerio', 'S.cerevisiae', 'E.coli', 'other'])
# species_group_order_reduced = np.array(['H.sapiens', 'M.musculus', 'D.melanogaster',
#                        'A.thaliana', 'S.cerevisiae', 'E.coli', 'rarely requested\norganisms'])
drop_species = ["D.rerio", "R.norvegicus"] 
rename_category = {'other': 'rarely requested\norganisms'}


def remove_species_and_rename_other(df, drop_species, rename_category):
    """
    Remove the species in drop_species, and rename the "other" category to 
    "rarely_requested_species".
    """
    reduced_df = df.loc[~df.species_group.isin(drop_species),:].copy()
    reduced_df.species_group.cat.remove_categories(drop_species, inplace = True)
    reduced_df.species_group.cat.rename_categories(rename_category, inplace = True)
    
    return reduced_df




def savefig_multiformat(filename, 
                        additional_formats=['png'], 
                        transparent = True,
                        dpi = 300,
                        *args, 
                        **kwargs):
    '''
    Save a figure in the file format specified by the file name
    as well as in a number of additional file formats.
    '''
    fig = plt.gcf()
    fig.savefig(filename, 
                transparent = transparent,
                dpi = dpi,
                *args,
                **kwargs)
    
    ## save in additional formats
    extension = filename.split('.')[-1]
    filename_wo_extension = filename.strip(extension)
    for ext in additional_formats:
        fig.savefig(filename_wo_extension+ext,
                    transparent = transparent,
                    dpi = dpi,
                    *args,
                    **kwargs)
        
        
        
        
def add_database_column(df, col, translate_to = 'database',
                        etype_to_database = etype_to_database,
                        database_to_etype = database_to_etype,
                        order = dbColors.index):
    '''
    In a dataframe that contains a column with etype or database, translate the column into the other and add it to the 
    dataframe. 
    
    order: The order of the database categories desired. None if no categorical conversion.
    '''

    original = df.loc[:,col].values

    if translate_to == 'database':
        database = [etype_to_database[e] for e in original]
        df = df.assign(database = database)

    elif translate_to == 'etype':
        etype = [database_to_etype[cat] for cat in original]
        df = df.assign(etype = etype)

    else:
        print("Choose 'database' or 'etype'.")
        return

    if order != None:
        df.loc[:,'database'] = pd.Categorical(df.database, categories = order, ordered = True)

    return df
