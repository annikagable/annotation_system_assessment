## For one species, plot the relationship between term sizes and the enrichment effect size (aka signal strength)

# In order to be able to plot the results from all annotation databases
# together in one plot, the PubMed terms were downsampled
# to 200'000 so that their size is comparable to the other
# annotation databases.

# # Warning:
# Although the files are called sigTermDf and dataId_isSig,
# The p-value threshold is here set to 1, which is why all
# terms with any overlap at all with the user input are 
# included here. 

## Usage: python scripts/plot_term_size_v_effect.py members_file termDf_file pval_plot_file  pval_plot_by_database_file effect_plot_file effect_plot_all_terms_file effect_plot_by_database_file sig_vs_insig_term_sizes_plot_file correlations_file

## Example: python scripts/plot_term_size_v_effect.py data/raw/global_enrichment_annotations/9606.terms_members.tsv data/results/cameraPR/downsampled_PubMed/aggregation/sigTermDf.tsv figures/cameraPR/9606.term_size_v_pval.svg  figures/cameraPR/9606.term_size_v_pval_by_database.svg figures/cameraPR/9606.term_size_v_effect_size.svg figures/cameraPR/9606.term_size_v_effect_size_all_terms.svg figures/cameraPR/9606.term_size_v_effect_size_by_database.svg figures/cameraPR/9606.sig_v_insig_term_size.svg data/results/term_size_v_effect_correlations.tsv


import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.ticker
import pandas as pd
import os
import seaborn as sns
import sys
import numpy as np
import warnings
import scipy.stats
import utils

def get_bin_means(data, x_col = 'log_term_size', y_col = 'abs_effect_size', n_bins = 100):
    '''
    Get the means of y_col (of dataframe data), binned by x_col.
    Required for plotting a smooth trend line.
    '''
    
    data.sort_values(by = x_col, inplace = True)
    
    x = data.loc[:, x_col]
    y = data.loc[:, y_col]
    
    # Get the mean y value per x bin
    counts,   bins = np.histogram(x, bins=n_bins)                 # count number of entries per bin
    bin_sums, bins = np.histogram(x, bins=n_bins, weights=y)      # sum of entries per bin
    bin_means = np.zeros((bin_sums.shape))+np.nan                 # mean value per bin: initialize
    condition = counts > 0                                        # mean value per bin: prevent zero divide
    bin_means[condition] = bin_sums[condition]/counts[condition]  # mean value per bin: sum/count
    
    # For plotting, get the midpoints of the bins
    bin_midpoints = []
    for i in range(len(bins)-1):
        lo = bins[i]
        hi = bins[i+1]
        bin_midpoints.append(np.mean([lo, hi]))
    bin_midpoints = np.array(bin_midpoints)
    
    # Also for plotting, remove the bins which did not have any entries.
    # This will prevent the connecting line from missing because either 
    # start or end point is missing
    bin_means = bin_means[counts!=0]
    bin_midpoints = bin_midpoints[counts!=0]
    
    return bin_midpoints, bin_means
        



def annotate_correlation(data, x_col, y_col, **kws):

    series_x = data.loc[:, x_col]
    series_y = data.loc[:, y_col]
    rho = scipy.stats.spearmanr(series_x, series_y)[0]
    ax = plt.gca()
    ax.text(.15, .8, f"ρ: {rho:.2f}",   
            transform=ax.transAxes, fontsize='small', 
            bbox=dict(facecolor='w', alpha=0.5))
    
def plot_trend(data, x_col, y_col, n_bins, **kws):
    bin_midpoints, bin_means = get_bin_means(data, x_col, y_col, n_bins)
    ax = plt.gca()
    bin_midpoints_linear = 10**bin_midpoints
    ax.plot(bin_midpoints_linear, bin_means, "-", color = "darkred", linewidth=2)
    
    
if __name__ == "__main__":
    
#     members_file = "data/raw/global_enrichment_annotations/9606.terms_members.tsv"
#     termDf_file  = "data/results/cameraPR/overlap_3-200/downsampled_PubMed/aggregation/sigTermDf_alpha1.tsv"

#     pval_plot_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_pval.svg"
#     pval_plot_by_database_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_pval_by_database.svg"

#     effect_plot_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_effect_size.svg"
#     effect_plot_all_terms_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_effect_size_all_terms.svg"
#     effect_plot_by_database_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_effect_size_by_database.svg"

#     effect_plot_FDR_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_effect_size_FDR.svg"
#     effect_plot_FDR_by_database_file = "figures/cameraPR/overlap_3-200/9606.term_size_v_effect_size_FDR_by_database.svg"
#     sig_vs_insig_term_sizes_plot_file = "figures/cameraPR/overlap_3-200/9606.sig_v_insig_term_size.svg"

#     correlations_file = "data/results/cameraPR/overlap_3-200/term_size_v_effect_correlations.tsv"

    assert(len(sys.argv) == 12)

    _, members_file, termDf_file,\
    pval_plot_file, pval_plot_by_database_file,\
    effect_plot_file, effect_plot_all_terms_file, effect_plot_by_database_file,\
    effect_plot_FDR_file, effect_plot_FDR_by_database_file, sig_vs_insig_term_sizes_plot_file,\
    correlations_file = sys.argv


    alpha = 0.05
    dbColors = utils.dbColors
    sns.set_context("notebook")
    talk = sns.plotting_context("talk")

    taxId = int(os.path.basename(members_file).split('.')[0])

    ## Read term sizes
    members = pd.read_table(members_file, usecols = [0,1,2], names = ["term_shorthand", "etype", "term_size"])

    ## Read all terms tested
    termDf = pd.read_table(termDf_file)
    termDf = termDf.loc[termDf.taxId == taxId, : ]
    termDf.loc[:,'database'] = pd.Categorical(termDf.database, categories=dbColors.index, ordered=True)

    ## Add term size and other columns to the term dataframe
    termDf = termDf.merge(members, how = "inner", on = "term_shorthand")
    termDf = termDf.assign(neg_log_pval = -np.log10(termDf.p_value),
                           abs_effect_size = abs(termDf.effect_size),
                           log_term_size = np.log10(termDf.term_size),
                           log_overlap = np.log10(termDf.overlap))

    sigPTermDf = termDf.loc[termDf.p_value < alpha, : ].copy()
    sigQTermDf = termDf.loc[termDf.q_value < alpha, : ].copy()

    rank_correlations = []
    
    if len(termDf) == 0:
        warnings.warn("termDf is empty")
    elif len(sigPTermDf) == 0:
        warnings.warn("sigPTermDf is empty")
    elif len(sigQTermDf) == 0:
        warnings.warn("sigQTermDf is empty")


    #####################################################
    #plot significant terms' sizes vs nonsignificant terms' sizes
    #####################################################

    termDf = termDf.assign(isSigQval = termDf.q_value <= alpha)

    ## Asssess if the term size difference between significantly enriched
    ## terms and non-significantly enriched terms is significant
    # def test_difference_in_distribution(df):
    #     x = df.loc[df.q_value <= alpha, 'log_term_size']
    #     y = df.loc[df.q_value >  alpha, 'log_term_size']
    #     U, pval = scipy.stats.mannwhitneyu(x, y)
    #     return pd.Series(dict(U=U,pval=pval))

    # gr = termDf.groupby("database")
    # WMW_test_df = gr.apply(test_difference_in_distribution)
    # WMW_test_df
    ## 	U	pval
    ## database		
    ## Reactome	3.119123e+10	0.000000e+00
    ## KEGG	1.939587e+09	1.839698e-88
    ## UniProt	2.916467e+09	0.000000e+00
    ## GO_BP	1.917490e+11	0.000000e+00
    ## GO_CC	1.129303e+10	0.000000e+00
    ## GO_MF	7.419540e+09	0.000000e+00
    ## SMART	5.900134e+08	0.000000e+00
    ## Pfam	1.953994e+09	0.000000e+00
    ## InterPro	7.355862e+09	0.000000e+00
    ## STRINGclusters	4.386382e+10	0.000000e+00
    ## PubMed	9.244340e+10	0.000000e+00
    ## => all are highly significant


    ## Plot significant terms vs insignificant terms grouped boxplot

    fig = plt.figure(figsize = (6,4.5))
    ax = fig.gca()
    if len(termDf) > 0:
        sns.boxplot(data =termDf, y = "database", hue= 'isSigQval', x = 'term_size', 
                      orient = 'h', showfliers = True, notch= True, fliersize = 1,
                      palette = ['lightgrey', 'tab:green']).set(xscale = 'log');
    plt.legend(title = "significantly\nenriched",
               fontsize= 'small', 
               title_fontsize='small', 
               loc = 'lower right');
    ax.set_ylabel('');
    ax.set_xlabel('term size');
    ax.set_xlim(0.9, 1.5e4)
    ax.set_title("Enriched annotation terms (H.sapiens)");
    fig.tight_layout()
    utils.savefig_multiformat(sig_vs_insig_term_sizes_plot_file)



    ##################################################
    ## term size v pval
    ##################################################

    with talk:
        alpha = alpha
        plt.figure(figsize = (6,4.5));
        if len(termDf) > 0:
            plt.hexbin(data = termDf, 
                       y = 'neg_log_pval',
                       x = 'term_size',
                       mincnt = 1,
                       norm = matplotlib.colors.LogNorm(),
                       xscale = 'log',
                       extent = (0, 4.4, 0, 45)
                      );

            plt.colorbar().set_label("number of terms");
        plt.xlabel("term size");
        plt.ylabel("-log10( raw p-value )");
        plt.hlines(y     = -np.log10(alpha),
                   xmin  = 0,
                   xmax  = 2.5118e4,
                   color = "grey");

    plt.tight_layout()
    utils.savefig_multiformat(pval_plot_file)



    ## Out of interest, we can also plot the mean pval across term size.
    ## This also shows a peak at term size 100, but doesn't drop so much towards
    ## the end, and gets increasingly wriggly for bigger term sizes, presumably
    ## because there are just a few such terms
    # bin_midpoints, bin_means = get_bin_means(data = termDf,#.loc[termDf.p_value < alpha], 
    #                                          x_col = 'log_term_size', 
    #                                          y_col = 'neg_log_pval', 
    #                                          n_bins = 100)
    # bin_midpoints_linear = 10**bin_midpoints
    # plt.plot(bin_midpoints_linear, bin_means, "-",
    #          color = "darkred", linewidth=2);  
    # plt.xscale('log')

    ##################################################
    ## term size v pval, facetted by database
    ##################################################

    with talk:
        g = sns.FacetGrid(termDf, col = "database", col_wrap = 4, 
                          sharey=True, sharex = False, height=4, aspect=1);
        if len(termDf) > 0:
            g = g.map(plt.hexbin, 'term_size', 'neg_log_pval', 
                      mincnt=1, gridsize=40, linewidths=0, 
                      xscale = 'log',
                      extent = (0, 4.4, 0, 50),
                      norm = matplotlib.colors.LogNorm());

            # Nice common colorbar
            cax = g.fig.add_axes([0.8, alpha, 0.03, 0.25]) # x0, y0, width, height
            cb = plt.colorbar(cax=cax)
            cb.set_label("# of terms")

        g.set_xlabels('term size')
        g.set_ylabels('-log10( raw p-value )')
        g.set_titles(col_template="{col_name}")




    g.fig.tight_layout()
    utils.savefig_multiformat(pval_plot_by_database_file)


    ##################################################
    ## term size v effect size, plot all terms, also non-enriched
    ##################################################
    bin_midpoints, bin_means = get_bin_means(data = termDf, 
                                             x_col = 'log_term_size', 
                                             y_col = 'abs_effect_size', 
                                             n_bins = 100)
    bin_midpoints_linear = 10**bin_midpoints

    rho, pval = scipy.stats.spearmanr(termDf.abs_effect_size, termDf.log_term_size)


    alpha = alpha
    with talk:
        fig = plt.figure(figsize = (6,4.5));
        ax = fig.gca()
        if len(termDf) > 0:
            hb = ax.hexbin(data = termDf, 
                       y = 'abs_effect_size',
                       x = 'term_size',
                       mincnt = 1,
                       xscale = 'log',
                       extent = (0, 4.4, 0, 1),
                       norm = matplotlib.colors.LogNorm(),
                      );
            fig.colorbar(hb).set_label("number of terms");
        ax.set_xlabel("term size");
        ax.set_ylabel("enrichment effect size");

        # Plot binned means (advantage: smoothens independent of point density)
        ax.plot(bin_midpoints_linear, bin_means, "-",
                 color = "darkred", linewidth=2);  

    #     # Plot correlation
    #     ax.text(.15, .4, f"Spearman's ρ: {rho:.2f}", transform=ax.transAxes)

    fig.tight_layout()
    utils.savefig_multiformat(effect_plot_all_terms_file)

    rank_correlations.append((effect_plot_all_terms_file, rho, pval))
    del rho

    ##################################################
    ## term size v effect size, enriched terms only
    ##################################################

    bin_midpoints, bin_means = get_bin_means(data = sigPTermDf, 
                                             x_col = 'log_term_size', 
                                             y_col = 'abs_effect_size', 
                                             n_bins = 100)
    bin_midpoints_linear = 10**bin_midpoints

    rho, pval = scipy.stats.spearmanr(sigPTermDf.abs_effect_size, sigPTermDf.log_term_size)
    # R   = scipy.stats.pearsonr( sigPTermDf.abs_effect_size, sigPTermDf.log_term_size)[0]

    with talk:
        alpha = alpha
        fig = plt.figure(figsize = (6,4.5));
        ax = fig.gca()
        if len(sigPTermDf) > 0:
            hb = ax.hexbin(data = sigPTermDf, 
                       y = 'abs_effect_size',
                       x = 'term_size',
                       mincnt = 1,
                       xscale = 'log',
                       extent = (0, 4.4, 0, 1),
                       norm = matplotlib.colors.LogNorm(),
                      );
            fig.colorbar(hb).set_label("number of terms");
        ax.set_xlabel("term size");
        ax.set_ylabel("enrichment effect size");

        # Plot binned means (advantage: smoothens independent of point density)
        ax.plot(bin_midpoints_linear, bin_means, "-",
                 color = "darkred", linewidth=2);  

    #     # Plot correlation
    #     ax.text(.2, .75, f"Pearson's R:\n  {R:.2f}", transform=ax.transAxes)
    #     ax.text(.2, .6, f"Spearman's rho:\n  {rho:.2f}", transform=ax.transAxes)

    fig.tight_layout()
    utils.savefig_multiformat(effect_plot_file)

    rank_correlations.append((effect_plot_file, rho, pval))
    del rho


    ##################################################
    ## term size v effect size, enriched terms only, facetted by database
    ##################################################

    x_col = 'term_size'
    x_col_log = 'log_term_size'
    y_col = 'abs_effect_size'

    with talk:
        g = sns.FacetGrid(sigPTermDf, col = "database", col_wrap = 4, 
                          sharey=False, sharex=False, height=4, aspect=1);
        if len(sigPTermDf) > 0:
            g = g.map(plt.hexbin, x_col, y_col, 
                      mincnt=1, gridsize=40, linewidths=0, 
                      xscale = 'log',
                      extent = (0, 4.4, 0, 1),
                      norm = matplotlib.colors.LogNorm());

            g.map_dataframe(annotate_correlation, x_col=x_col_log, y_col=y_col)
            g.map_dataframe(plot_trend,           x_col=x_col_log, y_col=y_col, n_bins=100)

            cax = g.fig.add_axes([0.8, alpha, 0.03, 0.25]) # x0, y0, width, height

            cb = plt.colorbar(cax=cax)
            cb.set_label("# of enriched terms")

        g.set_xlabels('term size')
        g.set_ylabels("enrichment effect size")
        g.set_titles(col_template="{col_name}")



    plt.tight_layout()   
    utils.savefig_multiformat(effect_plot_by_database_file)


    ##################################################
    ## term size v effect size, enriched terms only (FDR-based)
    ##################################################

    bin_midpoints, bin_means = get_bin_means(data = sigQTermDf, 
                                             x_col = 'log_term_size', 
                                             y_col = 'abs_effect_size', 
                                             n_bins = 100)
    bin_midpoints_linear = 10**bin_midpoints

    rho, pval = scipy.stats.spearmanr(sigQTermDf.abs_effect_size, sigQTermDf.log_term_size)
    # R   = scipy.stats.pearsonr( sigQTermDf.abs_effect_size, sigQTermDf.log_term_size)[0]


    with talk:
        fig = plt.figure(figsize = (6,4.5));
        ax = fig.gca()
        
        if len(sigQTermDf) > 0:
            hb = ax.hexbin(data = sigQTermDf, 
                       y = 'abs_effect_size',
                       x = 'term_size',
                       mincnt = 1,
                       xscale = 'log',
                       extent = (0, 4.4, 0, 1),
                       norm = matplotlib.colors.LogNorm(),
                      );
            fig.colorbar(hb).set_label("# of enriched terms");
        ax.set_xlabel("term size");
        ax.set_ylabel("enrichment effect size");

        # Plot binned means (advantage: smoothens independent of point density)
        ax.plot(bin_midpoints_linear, bin_means, "-",
                 color = "darkred", linewidth=2);  

    #     # Plot correlation 
    #     ax.text(.2, .75, f"Pearson's R:\n  {R:.2f}", transform=ax.transAxes)
    #     ax.text(.2, .6, f"Spearman's rho:\n  {rho:.2f}", transform=ax.transAxes)

    fig.tight_layout()
    utils.savefig_multiformat(effect_plot_FDR_file)

    rank_correlations.append((effect_plot_FDR_file, rho, pval))
    del rho

    ##################################################
    ## term size v effect size, enriched terms only (FDR-based), facetted by database
    ##################################################

    x_col = 'term_size'
    x_col_log = 'log_term_size'
    y_col = 'abs_effect_size'

    with talk:
        g = sns.FacetGrid(sigQTermDf, col = "database", col_wrap = 4, 
                          sharey=False, sharex=False, height=4, aspect=1);
        if len(sigQTermDf) > 0:
            g = g.map(plt.hexbin, x_col, y_col, 
                      mincnt=1, gridsize=40, linewidths=0, 
                      xscale = 'log',
                      extent = (0, 4.4, 0, 1),
                      norm = matplotlib.colors.LogNorm());
            
            cax = g.fig.add_axes([0.8, alpha, 0.03, 0.25]) # x0, y0, width, height
            cb = plt.colorbar(cax=cax)
            cb.set_label("# of enriched terms")


        g.map_dataframe(annotate_correlation, x_col=x_col_log, y_col=y_col)
        g.map_dataframe(plot_trend,           x_col=x_col_log, y_col=y_col, n_bins=100)

        g.set_xlabels('term size')
        g.set_ylabels("enrichment effect size")
        g.set_titles(col_template="{col_name}")


    g.fig.tight_layout()
    utils.savefig_multiformat(effect_plot_FDR_by_database_file)


    # ##################################################
    # ## overlap v effect size, enriched terms only (FDR-based)
    # ##################################################

    # bin_midpoints, bin_means = get_bin_means(data = sigQTermDf, 
    #                                          x_col = 'log_overlap', 
    #                                          y_col = 'abs_effect_size', 
    #                                          n_bins = 100)
    # bin_midpoints_linear = 10**bin_midpoints

    # rho, pval = scipy.stats.spearmanr(sigQTermDf.abs_effect_size, sigQTermDf.log_overlap)
    # #R   = scipy.stats.pearsonr( sigQTermDf.abs_effect_size, sigQTermDf.log_overlap)[0]

    #  
    #     fig = plt.figure(figsize = (6,4.5));
    #     ax = fig.gca()
    #     hb = ax.hexbin(data = sigQTermDf, 
    #                y = 'abs_effect_size',
    #                x = 'overlap',
    #                mincnt = 1,
    #                xscale = 'log',
    #                extent = (0, 4.4, 0, 1),
    #                norm = matplotlib.colors.LogNorm(),
    #               );
    #     fig.colorbar(hb).set_label("# of enriched terms");
    #     ax.set_xlabel("overlap size between term and user input");
    #     ax.set_ylabel("enrichment effect size");

    #     # Plot binned means (advantage: smoothens independent of point density)
    #     ax.plot(bin_midpoints_linear, bin_means, "-",
    #              color = "darkred", linewidth=2);  

    # #     # Plot correlation
    # #     ax.text(.2, .75, f"Pearson's R:\n  {R:.2f}", transform=ax.transAxes)
    # #     ax.text(.2, .6, f"Spearman's rho:\n  {rho:.2f}", transform=ax.transAxes)

    # fig.tight_layout()
    # fig.savefig(effect_plot_FDR_overlap_file)

    # rank_correlations.append((effect_plot_FDR_file, rho, pval))
    # del rho

    ###################################################
    ## Gather all rank correlations into one dataframe
    ###################################################

    correlation_df = pd.DataFrame.from_records(data = rank_correlations, 
                                               columns = ['plot_file_name', 
                                                          'Spearman_rho',
                                                          'pval'])

    correlation_df.to_csv(correlations_file, sep = '\t', index = False)







    # # aYH9Zi4CJh2g: many nans in the middle, 1900 non-nan entries
    # # p104qWt3rN6b: only 301 non-zero entries
    # # DFwDoizXLnN3: only 47 non-0.999999988956706 entries
    # # 1WaiwY4aHesu: only 149 non-zero entries
    # # PpsokxIpn2qU: only 173 non-zero entries
    # # OYvjazje8aAm: only 36 non-one entries



