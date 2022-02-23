## cameraPR is a functional enrichment method that does not return the effect size, only direction and p-value
## Here, we are reading the results from cameraPR and adding an effect size. This is calculated as:
##
## (observed deviation from the mean) / (maximum possible deviation from the mean)  
##
## Since cameraPR is based on ranks, we also use ranks to calculate this:
##
## (mean_term_rank - mean_input_rank) / (1 - mean_input_rank)
##
## A rank of 1 is the most extreme rank, which is why it is used here.
## Note that for simplicity, I am using zero-based ranks in the script below, so the 1 is actually a 0.
##
## Running this script take just between a few miliseconds to a few seconds (for PubMed), so it is pretty
## inefficient (but convenient) to start up a separate process for each combination of dataId with database.


def calculate_signal_strength(enrich_file, user_file, effect_size_file):

    ## Read the user input
    user_input = pd.read_table(user_file, header = None, 
                               index_col = 0, names = ["shorthand", "input_value"])
    user_input.sort_values(by= "input_value", inplace = True)
    user_input = user_input.assign(input_rank = np.array(range(len(user_input))))

    ## Create a dictionary that maps the genes to their ranks in the user input 
    shorthand_to_rank = user_input.input_rank.to_dict()



    ## Read enrichment for certain user input and certain database
    enrichment = pd.read_table(enrich_file, dtype = {"overlap_genes": str})

    def convert_string_to_list(gene_string):
        list_of_strings = gene_string.split()
        list_of_ints = [int(s) for s in list_of_strings]
        return list_of_ints

    enrichment = enrichment.assign(gene_list = enrichment.overlap_genes.apply(convert_string_to_list))



    ## Calculate the effect sizes of all terms in the enrichment table
    mean_input_rank = user_input.input_rank.mean()

    ranks_of_genes = enrichment.gene_list.apply(lambda shorthand_list: [shorthand_to_rank[sh] for sh in shorthand_list])
    mean_rank = ranks_of_genes.apply(np.mean)
    effect_size = mean_rank.apply(lambda mean_rank: (mean_rank - mean_input_rank) / (0 - mean_input_rank))

    enrichment = enrichment.assign(ranks_of_genes = ranks_of_genes,
                                   mean_rank = mean_rank,
                                   effect_size = effect_size)



    ## Write the file that includes the enrichment results plus effect sizes (but not all intermediate measures, unless this takes super long to compute)
    out_df = enrichment.loc[:,['term_shorthand', 'overlap',         'direction', 
                               'p_value',        'q_value',         'effect_size',
                               'overlap_genes']]
    out_df.to_csv(effect_size_file, sep = "\t", index=False)
    
    
    
    
    
import pandas as pd
import numpy as np
import os
import sys

_, enrich_file, user_file, effect_size_file = sys.argv   

calculate_signal_strength(enrich_file, user_file, effect_size_file)

