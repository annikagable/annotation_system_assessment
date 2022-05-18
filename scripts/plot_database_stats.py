## Plot term size, nr of terms per gene and % genome coverage
## call in py38_mpl35 environment
# When using the py38_mpl35 environment for plotting, affinity designer does not correctly
# import svg files, see bug report here https://github.com/matplotlib/matplotlib/issues/20910
# Therefore, I'm using pdf as output for now. The matplotlib version 3.5 gives much nicer looking plots than 3.1.

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.transforms as trans
import utils

_, term_coverage_file, genome_coverage_file, percent_genome_coverage_file, species_members_file, database_stats_file = sys.argv

## input

# taxId=9606
# term_size_threshold = 250
# lo = 0
# hi = int(term_size_threshold)
# term_size_limits = (lo, hi)

# term_coverage_file = f"data/results/database_stats/{taxId}.term_coverage_per_gene_{lo}-{hi}.tsv"
# genome_coverage_file = f"data/results/database_stats/{taxId}.genome_coverage_{lo}-{hi}.tsv"
# percent_genome_coverage_file = f"data/results/database_stats/{taxId}.percent_genome_coverage_{lo}-{hi}.tsv"
# species_members_file = "data/raw/global_enrichment_annotations/9606.terms_members.tsv"

## output
# database_stats_file = f"figures/database_stats/{taxId}.database_stats_{lo}-{hi}.pdf"

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 9
plt.rcParams['lines.markersize'] = 1
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.grid.axis'] = 'x'
plt.rcParams['grid.alpha'] = 0.5
plt.rcParams['xtick.bottom'] = False
plt.rcParams['patch.edgecolor'] = 'black'

dbColors = utils.dbColors
etype_to_database = utils.etype_to_database
database_names = utils.dbNames_break



def add_database_column(df, col, etype_to_database):
    
    etypes = df.loc[:,col].values
    
    # make sure the etype is given as a negative integer
    negative_etypes = [-abs(e) for e in etypes]
    
    database = [etype_to_database[e] for e in negative_etypes]
    df = df.assign(database = database)

    return df


species_members = pd.read_table(species_members_file, usecols = [0,1,2], names = ["term_id", "etype", "term_size"])
species_members_db = add_database_column(df = species_members, col = "etype", etype_to_database = etype_to_database)


term_coverage_df = pd.read_table(term_coverage_file)
#genome_coverage_df = pd.read_table(genome_coverage_file, index_col = 0)
percent_genome_coverage_df = pd.read_table(percent_genome_coverage_file, index_col = 0)


###### Plot number of terms and term size distribution for all databases ##########

# fig, axs = plt.subplots(1, 2, sharey=True, 
#                         figsize = (6,6), 
#                         gridspec_kw={'width_ratios': [3, 4]}
#                        )
# # Remove vertical space between axes
# # fig.subplots_adjust(wspace=0, hspace =0)

# # Plot each graph
# sns.countplot(data = species_members_db, y = "database", 
#               palette = dbColors, order = dbColors.index, 
#               log = True, ax = axs[0]);
# axs[0].set_ylabel("");
# axs[0].set_xlabel("# terms");

# sns.boxenplot(data = species_members_db, x = 'term_size', y = 'database', orient = 'h', ax = axs[1], 
#               palette = dbColors, order = dbColors.index);
# axs[1].tick_params(left = False)
# axs[1].set_xscale("log")
# axs[1].set_ylabel("");
# axs[1].set_xlabel("term size");
# plt.tight_layout(h_pad=0, w_pad=0)



####### Plot term and genome coverage ############################################

# fig, axs = plt.subplots(1, 2, sharey=True, figsize = (7,3), gridspec_kw={'width_ratios': [3, 4]})
# # Remove vertical space between axes
# #fig.subplots_adjust(wspace=0)

# # Plot 1: distribution per gene
# g = sns.boxenplot(data = term_coverage_df, y = "database", x = "term_count",
#                  palette = dbColors, order = dbColors.index, ax = axs[0]);
# g.set_xscale('symlog')
# g.set_ylabel('')
# g.set_xlabel('# of terms per gene')
# g.set_xlim(-0.5,None)

# # Plot 2: coverage
# cover_colors = '#b7308b #eb7655 #f8e326 #787878'.split()
# genome_coverage_df.plot(kind = 'barh', stacked = True, color = cover_colors, ax = axs[1]);
# #plt.legend(title = "# terms per gene", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon = False);
# plt.legend(title = "# terms\nper gene", framealpha = 1);
# plt.ylabel('');
# plt.xlabel('protein-coding genes');
# plt.tick_params(left = False)
# plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(base=5000));
# plt.gca().invert_yaxis()

# plt.tight_layout(h_pad=0, w_pad=0)



################ Plot all together, including nr of terms #########################################

# fig, axs = plt.subplots(1, 4, sharey=True, 
#                         figsize = (10,5), 
#                         #gridspec_kw={'width_ratios': [3, 4]}
#                        )
# fig.subplots_adjust(wspace=0, hspace = 0)
# # number of terms
# sns.countplot(data = species_members_db, y = "database", 
#               palette = dbColors, order = dbColors.index, 
#               log = True, ax = axs[0]);
# axs[0].set_ylabel("");
# axs[0].set_xlabel("# terms");

# # term sizes
# sns.boxenplot(data = species_members_db, x = 'term_size', y = 'database', orient = 'h', ax = axs[1], 
#               palette = dbColors, order = dbColors.index);
# axs[1].tick_params(left = False)
# axs[1].set_xscale("log")
# axs[1].set_ylabel("");
# axs[1].set_xlabel("term size");

# # terms per gene
# g = sns.boxenplot(data = term_coverage_df, y = "database", x = "term_count", ax = axs[2],
#                   palette = dbColors, order = dbColors.index);
# axs[2].tick_params(left = False)
# axs[2].set_xscale('symlog')
# axs[2].set_ylabel('')
# axs[2].set_xlabel('# of terms per gene')
# axs[2].set_xlim(-0.5,None)

# # genome coverage
# cover_colors = '#b7308b #eb7655 #f8e326 #787878'.split()
# genome_coverage_df.plot(kind = 'barh', stacked = True, color = cover_colors, ax = axs[3]);
# #plt.legend(title = "# terms per gene", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon = False);
# plt.legend(title = "# terms\nper gene", framealpha = 1);
# plt.tick_params(left = False)
# plt.xlabel('protein-coding genes');
# plt.tick_params(left = False)
# plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(base=5000));
# plt.gca().invert_yaxis()


################## Plot all together, except number of terms ###############################
mm = 1/25.4  # millimeter in inches
fig, axs = plt.subplots(1, 3, sharey=True, 
                        figsize = (8*21*mm,7*21*mm)
#                         figsize = (8,7)
                       )
fig.subplots_adjust(wspace=0, hspace = 0)

# term sizes
sns.boxenplot(data = species_members_db, x = 'term_size', y = 'database', 
              orient = 'h', ax = axs[0], 
              palette = dbColors, order = dbColors.index,
              linewidth = 1);
axs[0].tick_params(left = False)
axs[0].set_xscale("log")
axs[0].set_ylabel("");
axs[0].set_xlabel("term size");
# force the xtick positions, otherwise matplotlib makes too few
axs[0].xaxis.set_major_locator(ticker.FixedLocator([1e0,1e1,1e2,1e3,1e4]))
axs[0].set_xlim(7e-1,5e4)

# terms per gene
g = sns.boxenplot(data = term_coverage_df, y = "database", x = "term_count", 
                  ax = axs[1],
                  palette = dbColors, order = dbColors.index, linewidth = 1);
axs[1].tick_params(left = False)
axs[1].set_xscale('symlog')
axs[1].set_ylabel('')
axs[1].set_xlabel('terms per gene')
axs[1].set_xlim(-0.5,None)
# trying to bottom-align the 0 with the other xtick labels - doesn't work =(
# xtick_locations = axs[1].get_xticks()
# axs[1].set_xticks(ticks = xtick_locations)
# xticklabels = axs[1].get_xticklabels()
# axs[1].set_xticklabels(labels = xticklabels, va= "bottom")
# xticklabels[1].set_ha("left") #modify individual tick lables


# genome coverage
cover_colors = '#b7308b #eb7655 #f8e326 #787878'.split()
percent_genome_coverage_df.plot(kind = 'barh', stacked = True, color = cover_colors, ax = axs[2]);
plt.tick_params(left = False)
plt.legend(title = "# of terms per gene", bbox_to_anchor=(0.5, 1.05), 
           loc='lower center', borderaxespad=0., frameon = False);
plt.xlabel('protein-coding genome covered (%)');
plt.xlim(0, 100)
plt.gca().invert_yaxis()

ytick_locations = axs[2].get_yticks();
axs[2].set_yticks(ytick_locations);
# ## aligning the y database names doesn't work
axs[2].set_yticklabels(database_names, ha = "left");

# Get a tight bounding box around the figure without having to call plt.tight_layout,
# because that would create a small space between the subfigures (despite setting it to zero!)
tight_bbox_raw_x0y0 = np.array(axs[0].get_tightbbox(fig.canvas.get_renderer()))
tight_bbox_raw_x1y1 = np.array(axs[2].get_tightbbox(fig.canvas.get_renderer()))
tight_bbox_raw = trans.Bbox([tight_bbox_raw_x0y0[0], tight_bbox_raw_x1y1[1]])
tight_bbox = trans.TransformedBbox(tight_bbox_raw, trans.Affine2D().scale(1./fig.dpi))
plt.savefig(database_stats_file, dpi = 300, transparent = True, bbox_inches=tight_bbox, 
            metadata = {'Title': os.path.basename(database_stats_file)})
