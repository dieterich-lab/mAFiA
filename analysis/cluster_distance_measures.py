import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import kstest
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def get_homo_distances(in_df):
    all_distances = []
    for this_refseqID in in_df['refseqID'].unique():
        sub_df = in_df[in_df['refseqID'] == this_refseqID]
        if len(sub_df)>1:
            all_distances.extend(pdist(sub_df.rel_location.values[:, np.newaxis]))
    return all_distances


def get_hetero_distances(in_df_1, in_df_2):
    common_refseqIDs = list(set(in_df_1['refseqID'].unique()).intersection(set(in_df_2['refseqID'].unique())))
    all_distances = []
    for this_refseqID in common_refseqIDs:
        sub_df_1 = in_df_1[in_df_1['refseqID'] == this_refseqID]
        sub_df_2 = in_df_2[in_df_2['refseqID'] == this_refseqID]
        all_distances.extend(
            np.abs(sub_df_1.rel_location.values[:, np.newaxis] - sub_df_2.rel_location.values[np.newaxis, :]).flatten()
        )
    return all_distances


ds = 'HFpEF'
conditions = ['ctrl_merged', 'HFpEF_merged']

# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']

cond_names = {x: x.rstrip('_merged') for x in conditions}
cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['blue', 'red'])}

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
mod_colors = {'m6A': 'red', 'psi': 'purple'}

threshold = '50.0'

metaplot_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR'
img_out = '/home/adrian/img_out/site_cluster'
os.makedirs(img_out, exist_ok=True)

cond_mod_dist_measure = {
    this_cond: {} for this_cond in conditions}
for this_cond in conditions:
    for this_mod in mods:
        this_dist_measure_file = os.path.join(metaplot_dir, f"{ds}_{this_cond}_{this_mod}_modRatio{threshold}.dist.measures.txt")
        cond_mod_dist_measure[this_cond][this_mod] = pd.read_csv(this_dist_measure_file, sep='\t')

hist_range = [0, 3]
hist_num_bins = 50
bin_edges = np.linspace(*hist_range, hist_num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

ymax = 0.175
plt.figure(figsize=(10, 10))
for cond_ind, this_cond in enumerate(conditions):

    norm_histograms = {}

    plt.subplot(2, 2, cond_ind+1)
    hetero_dist = get_hetero_distances(cond_mod_dist_measure[this_cond]['m6A'],
                                       cond_mod_dist_measure[this_cond]['psi'])
    hetero_hist, _ = np.histogram(hetero_dist, bins=bin_edges)
    norm_hetero_hist = hetero_hist / np.sum(hetero_hist)
    norm_histograms['hetero'] = norm_hetero_hist
    for this_mod in mods:
        this_cond_mod_homo_dist = get_homo_distances(cond_mod_dist_measure[this_cond][this_mod])
        this_cond_mod_hist, _ = np.histogram(this_cond_mod_homo_dist, bins=bin_edges)
        norm_hist = this_cond_mod_hist / np.sum(this_cond_mod_hist)
        norm_histograms[this_mod] = norm_hist
        plt.plot(bin_centers, norm_hist, c=mod_colors[this_mod], label=rf'${{{dict_mod_display[this_mod]}}}$')
    plt.plot(bin_centers, norm_hetero_hist, c='g',
             label=rf"${{{dict_mod_display['m6A']}}}$" + r'$\longleftrightarrow$' + rf"${{{dict_mod_display['psi']}}}$")
    plt.legend()
    plt.ylim([0, ymax])
    plt.xticks([0, 1.5, 3], ['0', '50', '100'])
    plt.yticks(np.arange(0, ymax, 0.05))
    plt.xlabel('Transcript distance (%)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title(cond_names[this_cond], fontsize=15)
    # plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)

    plt.subplot(2, 2, 2+cond_ind+1)
    for this_mod in mods:
        plt.plot(bin_centers, np.cumsum(norm_histograms[this_mod]), c=mod_colors[this_mod],
                 label=rf'${{{dict_mod_display[this_mod]}}}$')
    plt.plot(bin_centers, np.cumsum(norm_histograms['hetero']), c='g',
             label=rf"${{{dict_mod_display['m6A']}}}$" + r'$\longleftrightarrow$' + rf"${{{dict_mod_display['psi']}}}$")
    res_ks = kstest(get_homo_distances(cond_mod_dist_measure[this_cond]['m6A']),
                    get_homo_distances(cond_mod_dist_measure[this_cond]['psi']))
    x_max_ks = bin_centers[np.argmax(np.abs(np.cumsum(norm_histograms['m6A']) - np.cumsum(norm_histograms['psi'])))]
    plt.axvline(x=x_max_ks, c='gray', alpha=0.5)
    plt.text(x_max_ks+0.05, 0.05, f'ks_stat = {res_ks.statistic:.3f}\np_val = {res_ks.pvalue:.3E}')
    plt.legend()
    plt.ylim([0, 1.05])
    plt.xticks([0, 1.5, 3], ['0', '50', '100'])
    plt.yticks(np.arange(0, 1.05, 0.25))
    plt.xlabel('Transcript distance (%)', fontsize=12)
    plt.ylabel('CDF', fontsize=12)
plt.suptitle(ds, fontsize=20)
plt.savefig(os.path.join(img_out, f'clustering_{ds}_with_ks_test.png'), bbox_inches='tight')