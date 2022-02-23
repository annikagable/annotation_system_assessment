import pandas as pd
import os
import sys


def write_filtered_deduplicated_data(filtered_user_inputs_file, non_duplicate_dataId_files, out_dir, table_out_file):
    '''
    I want two things: 
    a) a dataframe with all filtered data and an extra column for duplicate/non-duplicate
    b) each deduplicated user output written to a separate file, as input to the subsequent enrichment functions.
    '''
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    filtered_user_inputs = pd.read_table(filtered_user_inputs_file)
    filtered_user_inputs = filtered_user_inputs.assign(isDuplicate = True)
    #print(filtered_user_inputs.head())

    for dedup_dataId_file in non_duplicate_dataId_files:
        print(f"Reading {dedup_dataId_file}.")
        ## Read in the deduplicated dataIds for each species
        with open(dedup_dataId_file, 'r') as dedup:
            DATAIDS = [line.strip() for line in dedup.readlines()]

        ## Parse the taxId from the file name
        taxId = os.path.basename(dedup_dataId_file).split('.')[0]
        
        for dataId in DATAIDS:

            ## Write out each dataId's data into a separate file
            dedup = filtered_user_inputs.loc[filtered_user_inputs.dataId == dataId, ['shorthand', 'value']]
            dedup_file = os.path.join(out_dir, f'{dataId}.{taxId}.input.tsv')
            dedup.to_csv(dedup_file, sep = '\t', index = False, header = False)

            ## Add non-duplicate info to dataframe containing all filtered inputs
            filtered_user_inputs.loc[filtered_user_inputs.dataId == dataId, 'isDuplicate'] = False

    ## Write the dataframe to file
    filtered_user_inputs.to_csv(table_out_file, sep = "\t", index = False)
    


# filtered_user_inputs_file = "data/interim/filtered_user_inputs.tsv"
# TAXIDS = [7227, 10116]
# non_duplicate_dataId_files = [f"data/interim/deduplicated_dataIds/{taxId}.tsv" for taxId in TAXIDS]
# out_dir =  'data/interim/filtered_deduplicated_user_inputs'
# table_out_file = "data/interim/filtered_deduplicated_user_inputs.tsv"

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'w')
sys.sterr = sys.stdout

filtered_user_inputs_file = snakemake.input[0]
non_duplicate_dataId_files = snakemake.input[1:]
table_out_file = snakemake.output[0]
out_dir =  snakemake.output[1]
write_filtered_deduplicated_data(filtered_user_inputs_file, non_duplicate_dataId_files, out_dir, table_out_file)
