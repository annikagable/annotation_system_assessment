# After generating the data for the cumulative gene union size by term size 
# for summary statistics over all H.sapiens user inputs, for all possible 
# term sizes and for individual user's curves, I'm here plotting three user 
# inputs plus the summary into one single plot.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
import sys
import utils

_, cum_gene_unions_file, summary_cum_gene_unions_file, enrichment_specificity_file, enrichment_specificity_legend_file = sys.argv

# # input files
# # computed on all possible term sizes => for plotting summary curves
# summary_cum_gene_unions_file = "data/results/cameraPR/overlap_3-Inf/unique_enriched_genes/9606/mean_median_percentiles_of_cum_gene_union_by_term_size.tsv"
# # computed only on existing term sizes => for plotting individual users curves
# cum_gene_unions_file= "data/results/cameraPR/overlap_3-Inf/unique_enriched_genes/9606/cumulative_gene_unions_by_term_size.tsv"

# # output file
# enrichment_specificity_file = "figures/cameraPR/overlap_3-Inf/unique_enriched_genes/9606/enrichment_specificity.svg"
# enrichment_specificity_file = "figures/cameraPR/overlap_3-Inf/unique_enriched_genes/9606/enrichment_specificity.pdf"
# enrichment_specificity_legend_file = "figures/cameraPR/overlap_3-Inf/unique_enriched_genes/9606/enrichment_specificity_legend.svg"


# read input files
summary_by_term_size_df = pd.read_table(summary_cum_gene_unions_file)
cum_gene_unions_by_term_size_df = pd.read_table(cum_gene_unions_file)

sum_no_zero = summary_by_term_size_df.loc[summary_by_term_size_df.median_gene_union_size !=0, :]


dbDashes = utils.dbDashes
dbColors = utils.dbColors

# The context for saving the figure needs to be notebook, while the 
# actual plotting is done with "with sns.plotting_context('talk'):.
# Otherwise, the axis gets too few ticks.
sns.set_context('notebook')

def plot_user_curve(dataId, 
                    ax, 
                    example_number,
                    xlim,
                    ylim,
                    dataframe = cum_gene_unions_by_term_size_df, 
                    dbDashes = dbDashes,
                    dbColors = dbColors,
                    linewidth = 3,
                    yscale = "linear",
                    ):
    
    df = dataframe.loc[dataframe.dataId == dataId,:]

    g = sns.lineplot(data = df, x = "term_size",
                 y = "cumulative_gene_union_size", hue = 'database',
                 palette = dict(dbColors), legend = None, lw=linewidth,  
                 style = "database", dashes = dbDashes, ax = ax);    

    ax.set_xscale('log');
    ax.set_yscale(yscale);
    ax.set_xlim(xlim);
    ax.set_ylim(ylim);
    ax.set_title(f'example user query #{example_number}')
    ax.set_ylabel('# of unique enriched genes\n(cumulative)')
    ax.set_xlabel('term size')
    
    # Make sure that tick labels also appear when plotting as subplots
    # with sharey = True, sharex = True
    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    ax.yaxis.set_tick_params(which='both', labelleft=True)
    
    return g
            
            
def plot_median_curve(ax,
                      xlim,
                      ylim,
                      summary_by_term_size_df = summary_by_term_size_df,
                      dbDashes = dbDashes,
                      dbColors = dbColors,
                      yscale = "linear",
                      remove_zeros = False):
    alpha = 0.2
    
    # get columns containing the lower and upper bounds of the shaded confidence interval
    percentile_band_columns = [c for c in summary_by_term_size_df.columns if c.startswith("perc_")]
    assert len(percentile_band_columns) == 2
    
    
    
    if remove_zeros:
        summary_df = summary_by_term_size_df.loc[summary_by_term_size_df.median_gene_union_size !=0, :]
        condition = (summary_by_term_size_df.loc[:,percentile_band_columns[0]] !=0) &\
                    (summary_by_term_size_df.loc[:,percentile_band_columns[1]] !=0)
        percentile_band_df = summary_by_term_size_df.loc[condition, :]
    else:
        summary_df         = summary_by_term_size_df
        percentile_band_df = summary_by_term_size_df
        
    if len(summary_by_term_size_df) > 0:
        
        

        # lower bound to shaded area
        sns.lineplot(data = percentile_band_df, 
                     x = 'term_size', y = percentile_band_columns[0], 
                     linewidth = 0,
                     hue = 'database', hue_order=dbColors.index.to_list(), 
                     palette = dbColors.to_list(), alpha = alpha, legend = None,
                     ax = ax);

        # upper bound to shaded area
        sns.lineplot(data = percentile_band_df, 
                     x = 'term_size', y = percentile_band_columns[1], 
                     linewidth = 0,
                     hue = 'database', hue_order=dbColors.index.to_list(), 
                     palette = dbColors.to_list(), alpha = alpha, legend = None,
                     ax = ax);

        # create shaded area
        for db in dbColors.index:
            df = percentile_band_df.loc[percentile_band_df.database == db, :]
            ax.fill_between(df.term_size, 
                             df.loc[:, percentile_band_columns[0]], 
                             df.loc[:, percentile_band_columns[1]],
                             alpha=alpha, color = dbColors[db], linewidth = 0)

        # plot dashed lines
        g = sns.lineplot(data = summary_df, 
                         x = 'term_size', y = 'median_gene_union_size', 
                         linewidth = 3,
                         hue = 'database', hue_order=dbColors.index.to_list(), 
                         palette = dbColors.to_list(),legend = None,
                         style = "database", dashes = dbDashes, style_order= dbColors.index.to_list(),
                         ax = ax)
    
    
    ax.set_xscale('log')
    ax.set_yscale(yscale);
    ax.set_title("all queries, median", fontdict={'fontweight':'bold'})
    ax.set_ylabel('# of unique enriched genes per\nuser input (cumulative)')
    ax.set_xlabel('term size')
    ax.set_xlim(xlim);
    ax.set_ylim(ylim)
    
    # Make sure that tick labels also appear when plotting as subplots
    # with sharey = True, sharex = True
    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    ax.yaxis.set_tick_params(which='both', labelleft=True)
    
    
    
    
    
### plot main figure ###
xlim = (1,2.5e4) 
#ylim = (-10,695) 
ylim = (-10,560) 

with sns.plotting_context('talk'):
    fig, axs = plt.subplots(2, 2, sharey=True,sharex=True, figsize = (9.5,9))

    plot_user_curve(dataId='2JcTy5jJdg7A', ax=axs[0][0], example_number=1, xlim=xlim, ylim=ylim);
    plot_user_curve(dataId='wDmGkKCK7XEv', ax=axs[0][1], example_number=2, xlim=xlim, ylim=ylim);
    plot_user_curve(dataId='zz26G5go7Urs', ax=axs[1][0], example_number=3, xlim=xlim, ylim=ylim);
    plot_median_curve(ax= axs[1][1], xlim=xlim, ylim=ylim);
    plt.tight_layout()
    
#utils.savefig_multiformat(filename = enrichment_specificity_file) #not working for some reason
fig.savefig(fname  = enrichment_specificity_file, 
            transparent = True,
            dpi = 300,
            bbox_inches = "tight")


### plot main figure, log ###
xlim = (1,2.5e4) 
ylim = (1,695) 
yscale = "log"

with sns.plotting_context('talk'):
    fig, axs = plt.subplots(2, 2, sharey=True,sharex=True, figsize = (9.5,9))

    plot_user_curve(dataId='2JcTy5jJdg7A', ax=axs[0][0], example_number=1, xlim=xlim, ylim=ylim, yscale = yscale);
    plot_user_curve(dataId='wDmGkKCK7XEv', ax=axs[0][1], example_number=2, xlim=xlim, ylim=ylim, yscale = yscale);
    plot_user_curve(dataId='zz26G5go7Urs', ax=axs[1][0], example_number=3, xlim=xlim, ylim=ylim, yscale = yscale);
    plot_median_curve(ax= axs[1][1], xlim=xlim, ylim=ylim, yscale = yscale, remove_zeros = True);
    plt.tight_layout()
    
#utils.savefig_multiformat(filename = enrichment_specificity_file) #not working for some reason.
base, ext =  os.path.splitext(enrichment_specificity_file)
fig.savefig(fname  = f"{base}_{yscale}_no_zero_no_band_{ext}", 
            transparent = True,
            dpi = 300,
            bbox_inches = "tight")



### plot legend separately ###

labels = utils.dbNames

dbDashes_offset = {k:(0,dbDashes[k]) for k in dbDashes}
color_dashes = list(zip(dbColors.values, dbDashes_offset.values()))
custom_lines = [matplotlib.lines.Line2D([0], [0], 
                                        color=color, 
                                        lw=10, 
                                        linestyle = dash) for color, dash in color_dashes]
    

with sns.plotting_context('talk'):
    fig, ax = plt.subplots(figsize = (10,15))
    legend = ax.legend(custom_lines, labels, ncol = 1, handlelength = 10, 
                       labelspacing = 3, loc = 'center right',
                       frameon = False)
    ax.axis('off');
    for txt in legend.get_texts():
        txt.set_ha("left") # horizontal alignment of text item
        txt.set_x(-178) # x-position
        txt.set_y(20) # y-position


    plt.tight_layout()
utils.savefig_multiformat(filename = enrichment_specificity_legend_file)