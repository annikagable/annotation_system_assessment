import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys


def add_species_groups(df, n_grouped_species = 8):
    
    other_group_name = "rarely requested\n organisms"
    
    assert(all([col in df.columns for col in ["dataId", "species", "taxId"]]))
    
    ### Part 1: Determine the most common species
    
    ## Get each dataId and the associated species
    dataAndSpecies = df.loc[:,["dataId", "species", "taxId"]].drop_duplicates()
    
    ## Count how many dataIds per species we have
    species_counts = dataAndSpecies.groupby(['species', 'taxId']).dataId.count().sort_values(ascending = False)

    #import pdb; pdb.set_trace()
    
    ## Select the most frequent n species to be species groups
    species_group_long_names = species_counts.head(n_grouped_species).reset_index(level = 0).species#.drop(columns = 'dataId').species

    
    ### Part 2: Create the groups from the most common species
    ### This part could also be applied to any dataframe as long as we know the species to choose.
    
    ## Create shortened names using the initial of the genus plus the species name
    def create_short_name(species_name):
        name_list = species_name.split()
        initial = name_list[0][0]
        species = name_list[1]
        return f"{initial}.{species}"

    species_group_short_names = species_group_long_names.apply(create_short_name)
    species_group_short_names.name = "species_group"
    species_group_dict  = species_group_short_names.to_dict()

    ## If the species has no group, it's group is other
    def species_to_group(tax):
        return species_group_dict.get(tax, other_group_name)

    ## Add species group column
    df = df.assign(species_group = df.taxId.apply(species_to_group))
    
    ## Convert the species to categorical
    df.loc[:,'species_group'] = pd.Categorical(df.species_group, 
                                                    categories = species_group_short_names.to_list()+[other_group_name], 
                                                    ordered = True)
    
    return df



def generate_additional_file_formats(filename, formats = ['pdf', 'png']):
    '''
    For a given input file, return a list containing the original file name plus
    the file name with other extensions which are specified in formats.
    '''
    no_extension = os.path.splitext(filename)[0]
    additional_file_names = [no_extension+"."+ext for ext in formats]
    all_file_names = [filename] + additional_file_names
    return all_file_names
    
    

# ## input files
# user_inputs_file = "data/interim/filtered_deduplicated_user_inputs.tsv"
# species_taxIds_file = "data/raw/species.v11.0.txt"
# proteins_to_shorthands_file = "data/raw/proteins_to_shorthands.v11.tsv"
# hsa_protein_info_file = "data/raw/9606.protein.info.v11.0.txt.gz"
# proteome_sizes_file = "data/interim/taxId_to_proteome_size.tsv"

# script_name = "plot_user_input_statistics.py"

# ## output files
# input_groups_plot_file = "figures/input_analysis/input_count_and_input_size_by_species_group.svg"
# input_other_plot_file = "figures/input_analysis/input_count_and_input_size_for_other_species.svg"
# human_histogram_file = "figures/input_analysis/histogram_human_user_inputs_per_protein.svg"



assert(len(sys.argv) == 9)
script_name = sys.argv[0] 
input_files = sys.argv[1:6]
output_files = sys.argv[6:9]

user_inputs_file,\
species_taxIds_file,\
proteins_to_shorthands_file,\
hsa_protein_info_file,\
proteome_sizes_file = input_files

input_groups_plot_file, input_other_plot_file, human_histogram_file = output_files



talk = sns.plotting_context("talk")
whitegrid = sns.axes_style('whitegrid', 
                           rc= {'axes.edgecolor': '0',
                                'patch.force_edgecolor': False,
                                'xtick.bottom': True})

python_version = 'Python ' + sys.version.split()[0]
fig_metadata = {'Creator': python_version, 
                'Author': 'Annika Gable', 
                'Title': script_name}

## Russian palette from flatuicolors.com
Russian_palette = ["#ea8685", "#f5cd79", "#546de5", "#e15f41", "#c44569", "#786fa6", "#f78fb3", "#3dc1d3", "#596275", "#f19066"]


## read user inputs
user_inputs = pd.read_table(user_inputs_file)

## Select only the non-duplicate user inputs
user_inputs = user_inputs.loc[~user_inputs.isDuplicate,:]


## read translation table between species name and taxId
species_taxIds = pd.read_table(species_taxIds_file, 
                               usecols= [0,2], 
                               names = ["taxId", "species"], 
                               header=0)

## for histogram: read translation table between protein name and shorthand
proteins_to_shorthands = pd.read_table(proteins_to_shorthands_file,
                                      header = None, 
                                      names=['stringId', 'shorthand'])


## for histogram: read the full name and preferred name of each protein
hsa_protein_info = pd.read_table(hsa_protein_info_file, usecols = [0,1])
hsa_protein_info.columns = ["stringId", "preferred_name"]

## Get the proteome sizes
proteome_sizes = pd.read_table(proteome_sizes_file, index_col =0, squeeze=True)

## Add a species name column to the user inputs
user_inputs = user_inputs.merge(species_taxIds)

## Add a proteome size column
user_inputs = user_inputs.merge(proteome_sizes, left_on = 'taxId', right_index=True)

## Add a column for the species groups to the user inputs
user_inputs = add_species_groups(df = user_inputs, n_grouped_species = 8)

## Get the input sizes of each input
input_sizes_df = user_inputs.groupby(['species', 
                                      'species_group', 
                                      'dataId', 
                                      'proteome_size']).apply(len)
input_sizes_df = input_sizes_df.reset_index(name = 'input_size')

## Get the fraction of the query size / proteome size
input_size_fraction = input_sizes_df.input_size/input_sizes_df.proteome_size
input_sizes_df = input_sizes_df.assign(input_size_fraction = input_size_fraction)

######################################################################
## Plot number of queries and input sizes for all species groups
######################################################################


with talk:
    with whitegrid:
        fig, axs = plt.subplots(1, 2, sharey=True, figsize = (10,6))
        fig.subplots_adjust(wspace=0.1)

        order = input_sizes_df.species_group.cat.categories

        sns.countplot(data = input_sizes_df, 
                      y = "species_group", 
                      palette = Russian_palette, 
                      ax = axs[0], 
                      order = order)
        ## add number of user inputs
        alignment = 'right'
        for p in axs[0].patches:
            axs[0].annotate(s = f'{p.get_width()}', 
                            xy = (p.get_width()+2, p.get_y()+0.2), 
                            ha=alignment, va='top', color='black')
            alignment = 'left'
            
        axs[0].set_xscale('log')
        axs[0].set_xlim(10,2000)
        axs[0].set_ylabel("");
        axs[0].set_xlabel("# user queries");

        sns.boxenplot(data = input_sizes_df, 
                      y = 'species_group', 
                      x = 'input_size_fraction',
                      orient = 'h', 
                      palette = Russian_palette, 
                      showfliers = False,
                      ax = axs[1], 
                      order = order);
        #axs[1].axvline(2000, color = 'darkred', ls = '--')
        axs[1].set_ylabel("");
        axs[1].set_xlabel("query sizes / genome sizes");

output_files = generate_additional_file_formats(input_groups_plot_file)
for f in output_files:
    fig.savefig(f, metadata = fig_metadata, bbox_inches = 'tight', dpi = 300)    


######################################################################
## Plot number of queries and input sizes for "rarely requested organisms"
######################################################################


other_input_sizes = input_sizes_df.loc[input_sizes_df.species_group == "rarely requested\n organisms", :]


with talk:
    with whitegrid:
        fig, axs = plt.subplots(1, 2, sharey=True, figsize = (14,14))
        # Remove horizontal space between axes
        fig.subplots_adjust(wspace=0)

        order = other_input_sizes.species.value_counts().index
        
        # Plot only if there is data for rarely requested organisms
        if len(other_input_sizes) > 0:
            sns.countplot(data = other_input_sizes, y = "species", ax = axs[0], order = order)
            sns.boxplot(  data = other_input_sizes, y = 'species', ax = axs[1], order = order,
                          x = 'input_size_fraction', orient = 'h',  );
            
        axs[0].set_ylabel("");
        axs[0].set_xlabel("# user queries");
        axs[1].set_ylabel("");
        axs[1].set_xlabel("query sizes / genome sizes")


    fig.savefig(input_other_plot_file, metadata = fig_metadata, bbox_inches = 'tight', dpi = 300)



###################################################################
## Plot the query frequency histogram for human proteins (supplemental)
###################################################################

human_user_inputs = user_inputs.loc[user_inputs.taxId == 9606, :]
protein_counts = human_user_inputs.groupby("shorthand").apply(len)
protein_counts.name = "number_of_user_inputs"

# Get the HUGO IDs of each human shorthand by merging protein info of human
# with protein_to_shorthand of all species
protein_info_shorthands = proteins_to_shorthands.merge(hsa_protein_info, on = 'stringId')


# Add the HUGO names to the protein count table
named_protein_counts = protein_info_shorthands.merge(protein_counts, 
                                                     on = 'shorthand', 
                                                     how = "outer")
assert(len(named_protein_counts) == len(hsa_protein_info))

# Set proteins not found in the user inputs to zero
named_protein_counts.loc[named_protein_counts.number_of_user_inputs.isna(), 'number_of_user_inputs'] = 0
named_protein_counts.loc[:,"number_of_user_inputs"] = named_protein_counts.number_of_user_inputs.astype(int)

# Sort by how often a protein was in the user inputs
named_protein_counts.sort_values('number_of_user_inputs', inplace = True)

assert(all(named_protein_counts.columns == ['stringId', 'shorthand', 'preferred_name', 'number_of_user_inputs']))

with talk:
    plt.figure(figsize=(10, 5))
    g = sns.distplot(named_protein_counts.number_of_user_inputs, 
                     kde=False,
                     bins= 39,                          # Fixing the bins and the range allows
                     hist_kws={"range": (0, 780),       # me to have bins of width 20.
                               "rwidth":0.9});
    #g.set_title("Histogram of user input count per protein.");
    g.set_ylabel("count")
    g.set_xlabel("Number of user queries per gene (human only)")

    # Annotate the 3 most frequent genes
    for i in [-1,-2,-3]:
        preferred_name, nr_user_inputs = named_protein_counts.iloc[i, [2,3]].values
        g.annotate(preferred_name, 
                     xy=(nr_user_inputs, 100), 
                     xytext=(nr_user_inputs - abs(i) * 50, 800 - abs(i)*150), 
                     fontsize=15, 
                     arrowprops={'width':0.4,'headwidth':5,'color':'#333333'})

g.figure.savefig(human_histogram_file, metadata = fig_metadata, bbox_inches = 'tight', dpi = 300)
