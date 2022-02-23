# run in py38_mpl35 environment in order to get x_jitter

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import utils

def annotate_robust_correlation(data, x_col, y_col, **kws):
    rhos = []
    n= 10
    random_noises = np.random.uniform(low=-0.1, high=0.1, size=(n, len(data)))
    for i in range(n):
        series_x = data.loc[:, "term_count"] + random_noises[i]
        series_y = data.loc[:, "number_of_user_inputs"]
        rho, pval = scipy.stats.spearmanr(series_x, series_y)
        rhos.append(rho)
    mean_rho = np.mean(rhos)
    ax = plt.gca()
    ax.text(.75, .90, f"ρ: {mean_rho:.2f}",   
            transform=ax.transAxes, 
            #fontsize='x-large', bbox=dict(facecolor='w', alpha=1)
           )
    

result_dir = "data/results"
terms_per_gene_file  = os.path.join(result_dir, "database_stats/9606.term_coverage_per_gene_0-250.tsv")
inputs_per_gene_file = os.path.join(result_dir, "input_analysis/input_signal_median_imputed.tsv")

outfile = "figures/user_input_vs_database_coverage/interest_vs_coverage.svg"

# Read data (is human only)
terms_per_gene = pd.read_table(terms_per_gene_file)
inputs_per_gene = pd.read_table(inputs_per_gene_file)

# Merge the relevant columns of the two tables
tpg = terms_per_gene.loc[:,['term_count', 'stringId', 'preferred_name',
                            'protein_size', 'shorthand','database']].copy()
ipg = inputs_per_gene.loc[:,['shorthand', 'number_of_user_inputs']].copy()
terms_and_inputs_per_gene = tpg.merge(ipg)


# Plot
with sns.plotting_context("talk"):
    g = sns.lmplot(x = "term_count", y = "number_of_user_inputs", data = terms_and_inputs_per_gene, col = "database",
           col_wrap = 3, col_order = utils.dbColors.index,
           x_jitter = 0.2, fit_reg=False,
           facet_kws = {"sharex": False}, 
           scatter_kws = dict(alpha = 0.5, linewidth=0, s = 3),
           line_kws= dict(color="black")
            )
    g.set(xscale = 'symlog', xlim = (-1,None), #ylim = (None, 1188), 
          xlabel = "# of terms per gene", ylabel = "# of user inputs per gene"); 
    g.map_dataframe(annotate_robust_correlation, x_col="term_count", y_col="number_of_user_inputs");
    g.set_titles(col_template = '{col_name}')
    
utils.savefig_multiformat(outfile)
