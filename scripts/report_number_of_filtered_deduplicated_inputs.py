## Count the number of original user inputs, the ones which remain after filtering, and the
## ones remaining after deduplication, also by taxId

## Usage: report_number_of_filtered_deduplicated_inputs.py file_of_concatenated_user_inputs_infile submission_counts_outfile
##
## Example: report_number_of_filtered_deduplicated_inputs.py data/interim/filtered_deduplicated_user_inputs.tsv deduplicated_user_submission_counts_by_taxId.tsv

import pandas as pd
import sys

# Read the file where all filtered user submissions are there with their genes, values, dataIds, and whether the submission has
# been marked a duplicate (or near duplicate) or not.
if len(sys.argv) == 1:
    # we assume that we're in an interactive session
    filtered_user_inputs_file = "data/interim/filtered_deduplicated_user_inputs.tsv"
    submission_counts_by_taxId_outfile = "data/results/deduplicated_user_submission_counts_by_taxId.tsv"
else:    
    # we assume script execution on the command line
    filtered_user_inputs_file = sys.argv[1]
    submission_counts_by_taxId_outfile = sys.argv[2]

print(f"Reading {filtered_user_inputs_file} to determine number of filtered and near-duplicate inputs.")
filtered_user_inputs_df = pd.read_table(filtered_user_inputs_file)

# Keep only the dataId and the duplicate status of each user submission 
filtered_dataIds_df = filtered_user_inputs_df.loc[:, ["dataId", 'taxId', 'isDuplicate']].drop_duplicates()

# Count how many user submissions are in the table of filtered user submissions
count_after_filtering = len(filtered_dataIds_df)
print(f"Number of user submissions after filtering out small and infinite inputs: {count_after_filtering}")

# Count how many user submissions are not marked as duplicates
deduplicated_dataIds_df = filtered_dataIds_df.loc[filtered_dataIds_df.isDuplicate == False,:]
count_after_deduplication = len(deduplicated_dataIds_df)
print(f"Number of user submissions after removing highly correlated data: {count_after_deduplication}")

# Count how many non-duplicate user submissions there are by species
deduplicated_dataIds_series = deduplicated_dataIds_df.groupby('taxId').apply(len).sort_values(ascending = False)
deduplicated_dataIds_df = deduplicated_dataIds_series.to_frame(name = 'user_submission_count')

deduplicated_dataIds_df.to_csv(submission_counts_by_taxId_outfile, index = True, header = True, sep = "\t")
