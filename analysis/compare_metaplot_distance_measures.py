import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']

cond_names = {x: x.rstrip('_merged') for x in conditions}
cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['blue', 'red'])}

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresholds = ['0.0', '50.0']

metaplot_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR'
img_out = '/home/adrian/img_out/metatranscript'

# coverage_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/HFpEF'
# cond_coverage_files = {
#     'ctrl_merged': os.path.join(coverage_dir, 'ctrl_merged.geneBodyCoverage.txt'),
#     'HFpEF_merged': os.path.join(coverage_dir, 'ctrl_merged.geneBodyCoverage.txt'),
# }

# cond_coverage = {}
# for this_cond in conditions:
#     this_cond_df_coverage = pd.read_csv(cond_coverage_files[this_cond], sep='\t', header=None).T
#     this_cond_df_coverage = this_cond_df_coverage[1:]
#     this_cond_df_coverage = this_cond_df_coverage.rename(columns={0: 'percentile', 1: 'coverage'})
#     cond_coverage[this_cond] = np.float64(this_cond_df_coverage['coverage'].values)

cond_mod_thresh_dist_measure = {
    this_cond: {
        this_mod: {} for this_mod in mods
    } for this_cond in conditions}
for this_cond in conditions:
    for this_mod in mods:
        for this_thresh in thresholds:
            this_dist_measure_file = os.path.join(metaplot_dir, f"{ds}_{this_cond}_{this_mod}_modRatio{this_thresh}.dist.measures.txt")
            cond_mod_thresh_dist_measure[this_cond][this_mod][this_thresh] = pd.read_csv(this_dist_measure_file, sep='\t')

hist_range = [0, 3]
hist_num_bins = 50
bin_edges = np.linspace(*hist_range, hist_num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])


def smooth(y, box_pts=5):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    for this_cond in conditions:
        this_cond_mod_hist, _ = np.histogram(
            cond_mod_thresh_dist_measure[this_cond][this_mod][thresholds[1]].rel_location.values, bins=bin_edges
        )
        this_cond_mod_hist_all, _ = np.histogram(
            cond_mod_thresh_dist_measure[this_cond][this_mod]['0.0'].rel_location.values, bins=bin_edges
        )
        # norm_hist = this_cond_mod_hist / cond_coverage[this_cond]
        norm_hist = this_cond_mod_hist / this_cond_mod_hist_all
        norm_hist_smoothed = smooth(norm_hist)
        plt.plot(bin_centers, norm_hist_smoothed, c=cond_colors[this_cond], label=cond_names[this_cond])
    plt.legend()
    plt.axvline(x=1, c='gray', alpha=0.5)
    plt.axvline(x=2, c='gray', alpha=0.5)
    plt.xticks([0.5, 1.5, 2.5], ['5\' UTR', 'CDS', '3\' UTR'])
    plt.ylabel('$N_{S\geq50}$ / $N_{covered}$', fontsize=12)
    plt.ylim([0, 0.1])
    plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
plt.suptitle(ds, fontsize=20)
plt.savefig(os.path.join(img_out, f'profile_{ds}.png'), bbox_inches='tight')