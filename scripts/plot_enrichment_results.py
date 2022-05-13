## This script plots the two basic enrichment result metrics.
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import sys
import utils

# out_dir = "figures/cameraPR"
# sigTerm_file = "data/results/cameraPR/aggregation/sigTermDf.tsv"
# dataId_isSig_file = "data/results/cameraPR/aggregation/dataId_isSig.tsv"
# script_name = "plot_metrics.py"
# species_subset = "reduced"

script_name, sigTerm_file, dataId_isSig_file, out_dir, species_subset = sys.argv

python_version = 'Python ' + sys.version.split()[0]
fig_metadata = {'Creator': python_version, 'Author': 'Annika Gable', 'Title': script_name}

def as_png(filename):
    return filename.replace('.svg', '.png')

talk = sns.plotting_context('talk')
poster = sns.plotting_context('poster')

dbColors = utils.dbColors
drop_species = ["D.rerio", "R.norvegicus"]

# species_order = ["H.sapiens", "M.musculus", "R.norvegicus",
#                  "D.rerio", "D.melanogaster", "A.thaliana", 
#                  "S.cerevisiae", "E.coli", "other"]


## -----------------------------------------------------      
## Read and format enrichment data 


sigTermNamedGroupDf = pd.read_table(sigTerm_file)
dataId_isSig = pd.read_table(dataId_isSig_file)

other_organisms = "other"

if species_subset == "reduced":
    ## Remove some species
    sigTermNamedGroupDf = sigTermNamedGroupDf.loc[~sigTermNamedGroupDf.species_group.isin(drop_species),:]
    dataId_isSig = dataId_isSig.loc[~dataId_isSig.species_group.isin(drop_species),:]

    ## Rename "other" category
    other_organisms = 'rarely requested\norganisms'
    sigTermNamedGroupDf.loc[sigTermNamedGroupDf.species_group == 'other', 'species_group'] = other_organisms
    dataId_isSig.loc[dataId_isSig.species_group == 'other', 'species_group'] = other_organisms

## Get the total number of inputs (whether they got enrichment or not, per species_group)
input_counts = dataId_isSig.groupby("species_group").dataId.unique().apply(len)
input_counts.name = "all_input_counts"

# Re-shuffle input counts so that the "other" species are at the end
input_counts = input_counts.sort_values(ascending = False)
cond = (input_counts.index == other_organisms)
species = input_counts.loc[~cond]
other = input_counts.loc[cond]
input_counts = species.append(other)


# Convert database and species group columns to categorical so that the databases will be plotted in the right order
sigTermNamedGroupDf.loc[:, "database"] = pd.Categorical(sigTermNamedGroupDf.database,
                                                        categories = dbColors.index.to_list(), ordered = True) 
dataId_isSig.loc[:, "database"] = pd.Categorical(dataId_isSig.database,
                                                categories = dbColors.index.to_list(), ordered = True) 

sigTermNamedGroupDf.loc[:, "species_group"] = pd.Categorical(sigTermNamedGroupDf.species_group,
                                                             categories = input_counts.index.to_list(),
                                                             ordered = True) 
dataId_isSig.loc[:, "species_group"] = pd.Categorical(dataId_isSig.species_group,
                                                      categories = input_counts.index.to_list(),
                                                      ordered = True) 




## -----------------------------------------------------        
## How many significant terms per user input?

def count_sig_terms_per_user_input(sigTermNamedGroupDf, dbColors):
    important_columns = ['dataId', 'taxId', 'database', 'species_group']

    ## a dataframe containing the significant term counts for each database and each dataId
    countQByDataId = sigTermNamedGroupDf.groupby(important_columns,
                                                 as_index = False,
                                                 observed = True).q_value.count()
    countQByDataId.columns = important_columns + ['nr_sig_terms_per_input']

    ## add counts for the database-dataId combinations that did not get any significant terms

    # get all dataIds with their taxIds
    dataId_to_taxId = countQByDataId.loc[:, ['dataId', 'taxId', 'species_group']].drop_duplicates()

    # For each dataId (with taxId) and each database, create a row 
    databases = dbColors.index   
    combination_df_list = [dataId_to_taxId.assign(database= db) for db in databases]
    combination_df = pd.concat(combination_df_list, ignore_index = True)

    # Make sure that the df where sig terms are counted contains a row for each combination
    countQByDataId_all_combinations = countQByDataId.merge(combination_df, how = "outer")

    # Set the previously non-existing combinations to zero instead of nan
    condition = countQByDataId_all_combinations.nr_sig_terms_per_input.isna()
    countQByDataId_all_combinations.loc[condition, 'nr_sig_terms_per_input'] = 0

    countQByDataId_all_combinations.loc[:, "database"] = pd.Categorical(
                                                            countQByDataId_all_combinations.database,
                                                            categories = dbColors.index.to_list(), ordered = True) 
    return countQByDataId_all_combinations

countQByDataId_all_combinations = count_sig_terms_per_user_input(sigTermNamedGroupDf, dbColors)

    
    
with talk:
    g = plt.figure(figsize = (20,10))
    g = sns.boxplot(data = countQByDataId_all_combinations,
                    x = 'species_group',
                    y = 'nr_sig_terms_per_input',
                    hue = "database",
                    palette = dbColors, 
                    order = input_counts.index);
    plt.yscale('symlog')
    g.set_xlabel('')
    g.set_ylabel("Enrichment results – \n# of terms detected as significant per user input") 
    handles, labels = g.axes.get_legend_handles_labels()
    g.legend_.remove()
    g.legend(handles, labels, ncol=4, loc='upper right', frameon=True)
    
    ## Add the number of total user inputs
    y = -0.12
    g.text(x = -0.12, y = y, s= "# inputs", transform=g.transAxes)
    for i, short_name in enumerate(input_counts.index):
        formatted_count = f"{input_counts[short_name]:4}"
        g.text(x = 0.05+(1/len(input_counts))*i, y = y, s= formatted_count, transform=g.transAxes) 


figure_path = os.path.join(out_dir, 'nr_sig_terms_per_user_input.svg')
g.figure.savefig(figure_path, metadata = fig_metadata, bbox_inches = 'tight')
g.figure.savefig(as_png(figure_path), metadata = fig_metadata, bbox_inches = 'tight')
    
    
## -------------------------------------------------------------
## How many users get any kind of enrichment, by database, starting from the dataId_isSig df? 

sig_dataIds = dataId_isSig.loc[dataId_isSig.atLeastOneSigTerm, :]

with talk: #and sns.axes_style('whitegrid'):
    h = sns.catplot(data = sig_dataIds,
                kind = "count",
                y = "database",
                col = "species_group",
                palette = dbColors,
                sharex=False,
                height = 5,
                aspect = 0.4);
    
    h.set_titles(col_template = '{col_name}')
    
    for i, ax in enumerate(h.axes[0]):
        ax.axvline(input_counts[i], ls='--', label="total # of user inputs")
        ax.set_xlabel('')
        
#         #rotate label ## not working properly
#         labels = ax.get_xticklabels()
#         ax.set_xticklabels(labels, rotation=90, horizontalalignment='center');

    h.axes[0,0].set_ylabel('')

    h.fig.text(x = 0.5,
               y = 0,
               s= '# of user inputs with at least one significant term',
               horizontalalignment = "center")
    
figure_path = os.path.join(out_dir, 'at_least_one_significant_facetGrid.svg')
h.fig.savefig(figure_path, metadata = fig_metadata, bbox_inches = 'tight')
h.fig.savefig(as_png(figure_path), metadata = fig_metadata, bbox_inches = 'tight')
    
    

## The same as the above plot, but with switched x and y axes to be able to better see
## the performance difference between databases

with talk:
    species_groups = sig_dataIds.species_group.cat.categories.to_list()
    h = sns.catplot(data = sig_dataIds,
                kind = "count",
                x = "database",
                col = "species_group",
                palette = dbColors,
                sharey=False,
                height = 5,
                aspect = 0.4);
    
    h.set_titles('')#col_template = '{col_name}')
    
    
    for i, ax in enumerate(h.axes[0]):
        ax.axhline(input_counts[i], ls='--', label="total # of user inputs")
        ax.set_xlabel(species_groups[i], labelpad = 18)
        ax.set_xticklabels('')
        ax.set_xticks([])
        ax.yaxis.get_major_locator().set_params(integer=True)
        
    h.axes[0,0].set_ylabel('# of user inputs\nwith at least one significant term')

figure_path = os.path.join(out_dir, 'at_least_one_significant_facetGrid_vertical.svg')
h.fig.savefig(figure_path, metadata = fig_metadata, bbox_inches = 'tight')
h.fig.savefig(as_png(figure_path), metadata = fig_metadata, bbox_inches = 'tight')
    
    
    
# ## The data that the countplot procuces:

# g = sig_dataIds.groupby(['species_group', 'database']).apply(len)
# g.name = 'dataId_count'
# any_enrichment_count = g.reset_index()
