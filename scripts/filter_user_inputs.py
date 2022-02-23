## This script reads all protein-value user input data files in the specified input dir and puts them all into one big tsv.
## It adds columns for the unique input id and also the taxonomic identifier for each user entry in the tsv.
## This is not really utilizing snakemake's power in terms of parallelization, but this doesn't take very long.

import pandas as pd
import glob
import argparse
import os
import numpy as np

def _read_one_user_input(f, drop_undefined = True):

    df = pd.read_table(f, header = None, names = ["shorthand", "value"])
        
    if drop_undefined:
        # drop NaNs in case there are some in the input
        df.dropna(inplace = True)

        # drop infs in case there are some in the input
        df = df.loc[np.isfinite(df.value)]

    ## get unique identifier and organism id of the user input
    basename = os.path.basename(f)
    dataId, taxId, _, _, _ = basename.split('.')

    ## add dataId and organsm to dataframe
    df = df.assign(dataId = dataId, taxId = taxId)
    df.reset_index(inplace=True, drop=True)

    return df
    

def _get_ratio_of_most_frequent_values(df, top_n):
    '''
    Calculate the fraction top n most frequent values compared to the user input size.
    '''
    # How big is the user input?
    input_size = len(df)
    
    # How often does each value occur?
    value_frequencies = df.groupby('value').apply(len).sort_values()
    
    # How many unique input values are there?
    unique_value_count = len(value_frequencies)
    
    # Which are the n most frequent values?
    most_frequent_values = value_frequencies.tail(top_n).index.tolist()
    
    # How many genes have one of the n most frequent values?
    input_size_of_top_occurring_values = len(df.loc[ df.value.isin(most_frequent_values)])
    
    # What is the ratio of the n most frequent values to the total input size?
    ratio_most_frequent_values_to_input_size = input_size_of_top_occurring_values/input_size
    
    return ratio_most_frequent_values_to_input_size
    
    
    
    
    
def read_user_inputs(user_input_dir, search_for, min_finite_values, min_unique_values, max_ratio_top_value_to_input_size):
    '''
    Reading user inputs from file and merging them with the stringIds.
    
    user_input_dir: Directory of user input data
    min_finite_values: the minimum number of finite values a user input needs to have to be used. Currently using 500.
    min_unique_values: minimum number of uniquely occuring values a user input needs. Currently using 10.
    max_ratio_top_values_to_input_size: 
    search_for: the regex to get all the input files. e.g. "*.9606.input.normal.txt"
    
    Returns a list of user input frames
    '''
    
    # get user input
    files = glob.glob(os.path.join(user_input_dir, search_for))
    
    user_input_data = []
    
    for f in files:
        dataId = os.path.basename(f).split('.')[0]
        
        # Skip one specific input that weirdly enough submitted the string network file to 
        # global enrichment, which screws up effect size vs pval distribution.
        if dataId == 'Ibd31zo46YrA':
            print("Skipping dataId", dataId,  "because it is the string network file instead of actual data.")
            continue
            
        ## read file, dropping nans and infs
        df = _read_one_user_input(f)
        
        if len(df) < min_finite_values:
            print(f"Skipping dataId {dataId} because of less than {min_finite_values} finite values.")
            continue

        if len(df.value.unique()) < min_unique_values:
            print(f"Skipping dataId {dataId} because there are less than {min_unique_values} unique values.")
            continue
            
        # Skip cases with extremely large values like yScyQznwnhCV because then a mean 
        # cannot be computed by numpy within regular float type, would have to use long.
        if (df.value.max() > 1e200) or (df.value.min() < -1e200):
            print("Skipping dataId", dataId, "because of values > |1e200|.")
            continue
        
        # Check if not more than a certain percentage of the input consists of a single value
        ratio_most_frequent_value_to_input_size = _get_ratio_of_most_frequent_values(df, top_n = 1)
        if ratio_most_frequent_value_to_input_size > max_ratio_top_value_to_input_size:
            print(f"Skipping dataId {dataId} because more than {max_ratio_top_value_to_input_size*100}% of the input consists of a single value.")
            continue

        input_count = len(df)
                
        ## merge protein_info and input into one table => skipping for now, not necessary
        ## df = df.merge(protein_info)
        
        user_input_data.append(df)
    
    # Combine all user data into one dataframe
    user_input_df = pd.concat(user_input_data, ignore_index=True)
    
    return(user_input_df)

def write_out_separate(part_df):
    part_df.to_csv(f'{args.out_dir}/{part_df.taxId.values[0]}.{part_df.dataId.values[0]}',
    sep = '\t', 
    header = None, 
    index = False, 
    columns = ['shorthand', "value"])


### When calling this as a script, it will just read the user data all into one big dataframe with 
### a column for dataId and taxId added, and save to file.

if __name__ == "__main__":
    print("starting")
    parser=argparse.ArgumentParser()
    parser.add_argument('dir', help='input dir')
    parser.add_argument('--search_for', help='search string', default="*.input.normal.txt")
    parser.add_argument('--min_finite_values', help='How many valid entries does the user input have to have as a minimum requirement?', type=int, default = 500)
    parser.add_argument('--min_unique_values', help='How many unique entries does the user input have to have as a minimum requirement?', type=int, default = 10)
    parser.add_argument('--max_ratio_top_value_to_input_size', help='Which fraction of the input is allowed to consist of a single value?', type=float, default = 0.8)
    parser.add_argument('--out_file', help='output file name of the tsv containing all valid user inputs', default = "filtered_user_inputs.tsv")
    parser.add_argument('--out_dir', help='dir where the filtered user inputs are stored.', default = 'data/interim/filtered_user_inputs')    
    args=parser.parse_args()

    print(args)
    
    print("Reading user inputs.")
    user_input_df = read_user_inputs(user_input_dir = args.dir, 
                                     search_for = args.search_for,
                                     min_finite_values = args.min_finite_values,
                                     min_unique_values = args.min_unique_values,
                                     max_ratio_top_value_to_input_size = args.max_ratio_top_value_to_input_size)

    user_input_df.to_csv(args.out_file, sep = '\t', index = False)
    print(f"Saved to {args.out_file}")

#     os.mkdir(args.out_dir) 
#     user_input_df.groupby('dataId').apply(write_out_separate)
#     print(f"All filtered inputs saved to {args.out_dir}")

