import os
import pandas as pd
import numpy as np
import scipy.stats as stats

import matplotlib
matplotlib.use('PDF')
import mpl_toolkits.axisartist as axisartist
import matplotlib.pyplot as plt
import seaborn as sns

CHOSEN_REPLICATE = 57


def my_xformat(x):
    tx, ty = x.get_position()
    n = f"{tx:1.1e}"
    m, e = n.split('e')
    if e == '+00':
        t = f"${float(m):1.1f}$"
    elif '+' in e:
        t = f"${float(m):1.1f}\\times10^{{{e.lstrip('+0')}}}$"
    else:
        t = f"${float(m):1.1f}\\times10^{{-{e.lstrip('-0')}}}$"
    return matplotlib.text.Text(tx, ty, t)


def my_yformat(x):
    tx, ty = x.get_position()
    n = f"{ty:1.1e}"
    m, e = n.split('e')
    if e == '+00':
        t = f"${float(m):1.1f}$"
    elif '+' in e:
        t = f"${float(m):1.1f}\\times10^{{{e.lstrip('+0')}}}$"
    else:
        t = f"${float(m):1.1f}\\times10^{{-{e.lstrip('-0')}}}$"
    return matplotlib.text.Text(tx, ty, t)


if __name__ == '__main__':
    sns.set_style('ticks')

    palette = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a',
               '#66a61e', '#e6ab02', '#a6761d', '#666666']

    print('loading data...')

    d = pd.concat([pd.read_csv(f'../computed_matches/{f}')
                   for f in os.listdir('../computed_matches')])\
        .sort_values(['rep'])\
        .reset_index(drop=True)\
        .drop(columns=['s1', 's2'])

    """ FIGURE: MATCH DENSITIES """

    print('plotting match densities...')
    
    fig, axs = plt.subplots(1, 3, figsize=(
        12, 3), constrained_layout=True, gridspec_kw={'width_ratios': [27, 13, 20]})
    
    sns.histplot(d['true_allele_matches'], stat='density', element='poly', label='True',
                 discrete=True, fill=False, alpha=0.7, color=palette[0], ax=axs[0], marker='o', mfc='none', zorder=1)
    sns.histplot(d['called_allele_matches'], stat='density', element='poly',
                 label='Scheme 1\n(Called)', discrete=True, fill=False, alpha=0.7, color=palette[1], ax=axs[0], marker='+', mfc='none', zorder=2)
    sns.kdeplot(d['mean_allele_matches'],
                label='Scheme 2\n(Expected)', alpha=0.7, color=palette[2], ax=axs[0], zorder=0)
    
    axs[0].set_xlabel('Number of matches')
    axs[0].set_ylabel('Density')
    
    axs[0].xaxis.set_major_locator(plt.MultipleLocator(2))
    axs[0].xaxis.set_minor_locator(plt.MultipleLocator(1))
    axs[0].yaxis.set_major_locator(plt.MultipleLocator(0.02))
    
    sns.histplot(d['true_full_matches'], stat='density', element='poly', label='True',
                 discrete=True, fill=False, alpha=0.7, color=palette[0], ax=axs[1], marker='o', mfc='none', zorder=1)
    sns.histplot(d['called_full_matches'], stat='density', element='poly', label='Scheme 1\n(Called)',
                 discrete=True, fill=False, alpha=0.7, color=palette[1], ax=axs[1], marker='+', mfc='none', zorder=2)
    sns.kdeplot(d['mean_full_matches'], label='Scheme 2\n(Expected)',
                alpha=0.7, color=palette[2], ax=axs[1], zorder=0)
    
    axs[1].set_xlabel('Number of matches')
    axs[1].set_ylabel('Density')
    
    axs[1].xaxis.set_major_locator(plt.MultipleLocator(2))
    axs[1].xaxis.set_minor_locator(plt.MultipleLocator(1))
    
    sns.histplot(d['true_partial_matches'], stat='density', element='poly',
                 label='True', discrete=True, fill=False, alpha=0.7, color=palette[0], ax=axs[2], marker='o', mfc='none', zorder=1)
    sns.histplot(d['called_partial_matches'], stat='density', element='poly',
                 label='Scheme 1\n(Called)', discrete=True, fill=False, alpha=0.7, color=palette[1], ax=axs[2], marker='+', mfc='none', zorder=2)
    sns.kdeplot(d['mean_partial_matches'], label='Scheme 2\n(Expected)',
                alpha=0.7, color=palette[2], ax=axs[2], zorder=0)
    
    axs[2].set_xlabel('Number of matches')
    axs[2].set_ylabel('Density')
    
    axs[2].xaxis.set_major_locator(plt.MultipleLocator(2))
    axs[2].xaxis.set_minor_locator(plt.MultipleLocator(1))
    
    axs[2].legend()
    axs[2].get_legend().set_title('Matches')
    axs[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    axs[0].text(-0.15, 1.05, 'A', transform=axs[0].transAxes,
                size='large', fontweight='bold', va='baseline', ha='right')
    axs[1].text(-0.23, 1.05, 'B', transform=axs[1].transAxes,
                size='large', fontweight='bold', va='baseline', ha='right')
    axs[2].text(-0.19, 1.05, 'C', transform=axs[2].transAxes,
                size='large', fontweight='bold', va='baseline', ha='right')
    
    axs[0].text(0.5, 1.05, 'Allele matches', transform=axs[0].transAxes,
                size='large', fontweight='bold', va='baseline', ha='center')
    axs[1].text(0.5, 1.05, 'Fully matching loci', transform=axs[1].transAxes,
                size='large', fontweight='bold', va='baseline', ha='center')
    axs[2].text(0.5, 1.05, 'Partially matching loci', transform=axs[2].transAxes,
                size='large', fontweight='bold', va='baseline', ha='center')
    
    axs[0].set_xlim(-1, 26)
    axs[0].set_ylim(0, 0.181)
    axs[1].set_xlim(-1, 12)
    axs[2].set_xlim(-1, 19)
    
    fig.align_labels()
    
    fig.savefig('../figures/match_densities.pdf', bbox_inches='tight')

    """ FIGURE: MATCH REGRESSIONS """

    print('plotting match regressions...')

    fig, axs = plt.subplots(2, 3, figsize=(10, 6), constrained_layout=True)

    pal = sns.light_palette(palette[2], as_cmap=True)
    # pal = sns.cubehelix_palette(start=244*3/360, light=1, dark= 0.8, rot=0.4, gamma=3, as_cmap=True)
    # pal = sns.cubehelix_palette(start=1, rot=-0.3, gamma=1.2, as_cmap=True)

    #d_r = d[d.rep == CHOSEN_REPLICATE]
    d_r = d

    # figure out the range of values for the colorbar
    ca_hist, ca_x, ca_y = np.histogram2d(d_r['true_allele_matches'], d_r['called_allele_matches'], bins=25, range=[[0, 26], [0, 26]], density=False)
    cf_hist, cf_x, cf_y = np.histogram2d(d_r['true_full_matches'], d_r['called_full_matches'], bins=9, range=[[0, 12], [0, 12]], density=False)
    cp_hist, cp_x, cp_y = np.histogram2d(d_r['true_partial_matches'], d_r['called_partial_matches'], bins=19, range=[[0, 19], [0, 19]], density=False)
    ma_hist, ma_x, ma_y = np.histogram2d(d_r['true_allele_matches'], d_r['mean_allele_matches'], bins=25, range=[[0, 25], [0, 25]], density=False)
    mf_hist, mf_x, mf_y = np.histogram2d(d_r['true_full_matches'], d_r['mean_full_matches'], bins=9, range=[[0, 12], [0, 12]], density=False)
    mp_hist, mp_x, mp_y = np.histogram2d(d_r['true_partial_matches'], d_r['mean_partial_matches'], bins=19, range=[[0, 19], [0, 19]], density=False)

    ca_max = np.max(ca_hist)
    cf_max = np.max(cf_hist)
    cp_max = np.max(cp_hist)
    ma_max = np.max(ma_hist)
    mf_max = np.max(mf_hist)
    mp_max = np.max(mp_hist)

    ca_hist = np.ma.masked_equal(ca_hist, 0)
    cf_hist = np.ma.masked_equal(cf_hist, 0)
    cp_hist = np.ma.masked_equal(cp_hist, 0)
    ma_hist = np.ma.masked_equal(ma_hist, 0)
    mf_hist = np.ma.masked_equal(mf_hist, 0)
    mp_hist = np.ma.masked_equal(mp_hist, 0)

    # print(ca_max, cf_max, cp_max, ma_max, mf_max, mp_max)
    maxval = np.max([ca_max, cf_max, cp_max, ma_max, mf_max, mp_max])
    norm = matplotlib.colors.Normalize(vmin=1, vmax=3_000_000)
    # norm = matplotlib.colors.LogNorm(vmin=1, vmax=maxval)

    axs[0,0].pcolormesh(ca_x, ca_y, ca_hist.T, cmap=pal, norm=norm)
    axs[0,1].pcolormesh(cf_x, cf_y, cf_hist.T, cmap=pal, norm=norm)
    axs[0,2].pcolormesh(cp_x, cp_y, cp_hist.T, cmap=pal, norm=norm)
    axs[1,0].pcolormesh(ma_x, ma_y, ma_hist.T, cmap=pal, norm=norm)
    axs[1,1].pcolormesh(mf_x, mf_y, mf_hist.T, cmap=pal, norm=norm)
    axs[1,2].pcolormesh(mp_x, mp_y, mp_hist.T, cmap=pal, norm=norm)

    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=pal), ax=axs[:,-1], label='Number of pairs')
    cbar.set_ticks(np.arange(0, 3_000_001, 500_000))
    cbar.set_ticklabels(['0', '$5\\times 10^5$', '$1\\times 10^6$', '$1.5\\times 10^6$', '$2\\times 10^6$', '$2.5\\times 10^6$', '$\geq 3\\times 10^6$'])

    for ax in axs.flatten():
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.set_aspect('equal')

    for ax in axs[:, 0]:
        ax.set_xlim(0, 26)
        ax.set_ylim(0, 26)
        ax.xaxis.set_major_locator(plt.MultipleLocator(6))
        ax.yaxis.set_major_locator(plt.MultipleLocator(6))

        ax.plot([-1, 25], [-1, 25], color='black', linestyle='--')

    for ax in axs[:, 1]:
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 12)
        ax.xaxis.set_major_locator(plt.MultipleLocator(2))
        ax.yaxis.set_major_locator(plt.MultipleLocator(2))

        ax.plot([0, 9], [0, 9], color='black', linestyle='--')

    for ax in axs[:, 2]:
        ax.set_xlim(0, 19)
        ax.set_ylim(0, 19)
        ax.xaxis.set_major_locator(plt.MultipleLocator(2))
        ax.yaxis.set_major_locator(plt.MultipleLocator(2))

        ax.plot([0, 19], [0, 19], color='black', linestyle='--')

    for ax in axs[-1, :]:
        ax.set_xlabel('True matches')

    axs[0, 0].set_ylabel('Called matches')
    axs[1, 0].set_ylabel('Expected matches')

    # just show spearman correlation as $\rho$
    rho = stats.spearmanr(d_r.true_allele_matches,
                          d_r.called_allele_matches).correlation
    axs[0, 0].text(0.05, 0.9, '$\\rho$ = {:.2f}'.format(
        rho), transform=axs[0, 0].transAxes)

    rho = stats.spearmanr(d_r.true_allele_matches,
                          d_r.mean_allele_matches).correlation
    axs[1, 0].text(0.05, 0.9, '$\\rho$ = {:.2f}'.format(
        rho), transform=axs[1, 0].transAxes)

    rho = stats.spearmanr(d_r.true_full_matches,
                          d_r.called_full_matches).correlation
    axs[0, 1].text(0.05, 0.9, '$\\rho$ = {:.2f}'.format(
        rho), transform=axs[0, 1].transAxes)

    rho = stats.spearmanr(d_r.true_full_matches,
                          d_r.mean_full_matches).correlation
    axs[1, 1].text(0.05, 0.9, '$\\rho$ = {:.2f}'.format(
        rho), transform=axs[1, 1].transAxes)

    rho = stats.spearmanr(d_r.true_partial_matches,
                          d_r.called_partial_matches).correlation
    axs[0, 2].text(0.05, 0.9, '$\\rho$ = {:.2f}'.format(
        rho), transform=axs[0, 2].transAxes)

    rho = stats.spearmanr(d_r.true_partial_matches,
                          d_r.mean_partial_matches).correlation
    axs[1, 2].text(0.05, 0.9, '$\\rho$ = {:.2f}'.format(
        rho), transform=axs[1, 2].transAxes)

    axs[0, 0].text(-0.05, 1.05, 'A', transform=axs[0, 0].transAxes,
                   fontsize=12, fontweight='bold', va='baseline', ha='right')
    axs[0, 1].text(-0.05, 1.05, 'B', transform=axs[0, 1].transAxes,
                   fontsize=12, fontweight='bold', va='baseline', ha='right')
    axs[0, 2].text(-0.05, 1.05, 'C', transform=axs[0, 2].transAxes,
                   fontsize=12, fontweight='bold', va='baseline', ha='right')
    axs[1, 0].text(-0.05, 1.05, 'D', transform=axs[1, 0].transAxes,
                   fontsize=12, fontweight='bold', va='baseline', ha='right')
    axs[1, 1].text(-0.05, 1.05, 'E', transform=axs[1, 1].transAxes,
                   fontsize=12, fontweight='bold', va='baseline', ha='right')
    axs[1, 2].text(-0.05, 1.05, 'F', transform=axs[1, 2].transAxes,
                   fontsize=12, fontweight='bold', va='baseline', ha='right')

    axs[0, 0].text(0.5, 1.05, 'Allele matches', transform=axs[0, 0].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='center')
    axs[0, 1].text(0.5, 1.05, 'Fully matching loci', transform=axs[0, 1].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='center')
    axs[0, 2].text(0.5, 1.05, 'Partially matching loci', transform=axs[0,
                   2].transAxes, size='large', fontweight='bold', va='baseline', ha='center')
 
    axs[0, 0].text(-0.32, 0.5, 'Scheme 1', ha='center', va='center',
                   size='large', fontweight='bold', rotation=90, transform=axs[0, 0].transAxes)
    axs[1, 0].text(-0.32, 0.5, 'Scheme 2', ha='center', va='center',
                   size='large', fontweight='bold', rotation=90, transform=axs[1, 0].transAxes)

    fig.align_labels()

    # fig.savefig('figures/matches_regression.pdf', bbox_inches='tight')
    fig.savefig('../figures/matches_correlation.pdf',
                dpi=600, bbox_inches='tight')

    """ FIGURE: WILCOXON EFFECT SIZES """

    print('plotting effect sizes...')

    d_wilcoxon = pd.read_csv('wilcoxon_tests.csv')

    fig = plt.figure(figsize=(9, 5), constrained_layout=True)

    axs = np.zeros((2, 3), dtype=object)
    for i in range(2):
        for j in range(3):
            axs[i, j] = fig.add_subplot(
                2, 3, 3*i+j+1, axes_class=axisartist.Subplot)

    col = '#1b9e77'

    sns.histplot(d_wilcoxon.tc_allele_medians,
                 color=col, binwidth=0.00001, ax=axs[0, 0])
    sns.histplot(d_wilcoxon.tc_full_medians, color=col, binwidth=0.00001, ax=axs[0, 1])
    sns.histplot(d_wilcoxon.tc_partial_medians,
                 color=col, binwidth=0.00001, ax=axs[0, 2])

    sns.histplot(d_wilcoxon.tm_allele_medians,
                 color=col, binwidth=0.01, ax=axs[1, 0])
    sns.histplot(d_wilcoxon.tm_full_medians, color=col, binwidth=0.01, ax=axs[1, 1])
    sns.histplot(d_wilcoxon.tm_partial_medians,
                 color=col, binwidth=0.01, ax=axs[1, 2])

    for ax in axs.flatten():
        ax.set_xlabel('')
        ax.set_ylabel('')

    for ax in axs[1, :]:
        ax.set_xlabel('Estimated difference')

    for ax in axs[:, 0]:
        ax.set_ylabel('Count')

    axs[0, 0].text(0.5, 1.05, 'Allele matches', transform=axs[0, 0].transAxes,
                   ha='center', va='baseline', fontweight='bold', size='large')
    axs[0, 1].text(0.5, 1.05, 'Fully matching loci', transform=axs[0, 1].transAxes,
                   ha='center', va='baseline', fontweight='bold', size='large')
    axs[0, 2].text(0.5, 1.05, 'Partially matching loci', transform=axs[0,
                   2].transAxes, ha='center', va='baseline', fontweight='bold', size='large')

    # axs[0, 0].text(-0.3, 0.5, 'True vs. BEAGLE-called', transform=axs[0, 0].transAxes,
    #                size='large', ha='center', va='center', rotation=90, fontweight='bold')
    # axs[1, 0].text(-0.3, 0.5, 'True vs. BEAGLE-AP', transform=axs[1, 0].transAxes,
    #                size='large', ha='center', va='center', rotation=90, fontweight='bold')
    axs[0, 0].text(-0.3, 0.5, 'Scheme 1', ha='center', va='center',
                   size='large', fontweight='bold', rotation=90, transform=axs[0, 0].transAxes)
    axs[1, 0].text(-0.3, 0.5, 'Scheme 2', ha='center', va='center',
                   size='large', fontweight='bold', rotation=90, transform=axs[1, 0].transAxes)

    for ax in axs.flatten():
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(-max(abs(xmin), abs(xmax)), max(abs(xmin), abs(xmax)))

    # for ax in axs.flatten():
    #     ticks = ax.get_xticks().tolist()
    #     ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
    #     ax.set_xticklabels([my_xformat(t) for t in ax.get_xticklabels()])

    ## FORCE SAME Y SCALE IN ROWS, following 3rd plot in first row, first plot in second row
    max_xlim_0 = axs[0,2].get_xlim()[1]
    max_xlim_1 = axs[1,0].get_xlim()[1]
    max_xticks_0 = axs[0,2].get_xticks().tolist()
    max_xticks_1 = axs[1,0].get_xticks().tolist()

    for ax in axs[0,:]:
        ax.set_xlim(-max_xlim_0, max_xlim_0)
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(max_xticks_0))
        ax.set_xticklabels([my_xformat(t) for t in ax.get_xticklabels()])

    for ax in axs[1,:]:
        ax.set_xlim(-max_xlim_1, max_xlim_1)
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(max_xticks_1))
        # ax.set_xticklabels([my_xformat(t) for t in ax.get_xticklabels()])
        
    for ax in axs.flatten():
        ax.axis["bottom"].major_ticklabels.set(fontsize=9, va="baseline")
        ax.axis["right"].major_ticks.set_visible(False)
        ax.axis["top"].major_ticks.set_visible(False)

    # find max ylim and set it for all panels
    max_ylim = max([ax.get_ylim()[1] for ax in axs.flatten()])
    for ax in axs.flatten():
        ax.set_ylim(0, max_ylim) 


    axs[0, 0].text(-0.05, 1.05, 'A', transform=axs[0, 0].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='right')
    axs[0, 1].text(-0.05, 1.05, 'B', transform=axs[0, 1].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='right')
    axs[0, 2].text(-0.05, 1.05, 'C', transform=axs[0, 2].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='right')
    axs[1, 0].text(-0.05, 1.05, 'D', transform=axs[1, 0].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='right')
    axs[1, 1].text(-0.05, 1.05, 'E', transform=axs[1, 1].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='right')
    axs[1, 2].text(-0.05, 1.05, 'F', transform=axs[1, 2].transAxes,
                   size='large', fontweight='bold', va='baseline', ha='right')

    fig.align_labels()

    fig.savefig('../figures/wilcoxon_effect_sizes.pdf', bbox_inches='tight')

    print('done!')
