### Take input data which already has been cleaned for NA's and infinite values, as well as too large values,
### and remove inputs which are too similar in terms of correlation or the protein set studied. 

import pandas as pd
import numpy as np
import os
import sys

def pivot_to_species_matrix_list(all_data, out_dir):
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    df_dict = dict(list((all_data.groupby('taxId'))))
    mat_dict = dict()
    for taxId in df_dict:
        print(f"{taxId=}")
        data_matrix = df_dict[taxId].pivot_table(index = 'dataId', columns='shorthand', values='value', fill_value=np.nan)
        
        out_file = os.path.join(out_dir, f"{taxId}.tsv")
        data_matrix.to_csv(out_file, sep = '\t', index = True, header = True)
        
#         mat_dict[taxId] = data_matrix
#     return mat_dict

_, filtered_inputs, out_dir = sys.argv
print(f"{filtered_inputs =}") 
print(f"{out_dir =}")

## Read in filtered data
all_data = pd.read_table(filtered_inputs)

## Pivot long-form input data into a dictionary with one dataIds X shorthands dataframe per species.
pivot_to_species_matrix_list(all_data = all_data, out_dir = out_dir)