import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def smooth(y, box_pts=5):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


# ds = 'HFpEF'
# sub_ds = 'ctrl'
# sub_ds = 'HFpEF'

# ds = 'Diet'
# sub_ds = 'WT_CD'
# sub_ds = 'WT_WD'

# ds = 'TAC'
# sub_ds = 'SHAM'
# sub_ds = 'TAC'

ds = 'CM'
conditions = ['WT', 'M3KO']

cond_names = {this_cond: this_cond for this_cond in conditions}

# conditions = ['below_90', 'above_90']
# conditions = ['below_50', 'above_150']
# cond_names = {this_cond: this_cond.replace('_', ' ') + 'nts' for this_cond in conditions}

cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['blue', 'red'])}

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresholds = ['0.0', '50.0']

metaplot_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR'
img_out = '/home/adrian/img_out/metatranscript'

# metaplot_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA'
# img_out = '/home/adrian/img_out/metatranscript_polyA'
os.makedirs(img_out, exist_ok=True)

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
hist_num_bins = 60
bin_edges = np.linspace(*hist_range, hist_num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

plt.figure(figsize=(10, 4))
all_norm_hist = []
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    for this_cond in conditions:
        this_cond_mod_hist, _ = np.histogram(
            cond_mod_thresh_dist_measure[this_cond][this_mod][thresholds[1]].rel_location.values, bins=bin_edges
        )
        this_cond_mod_hist_all, _ = np.histogram(
            cond_mod_thresh_dist_measure[this_cond][this_mod]['0.0'].rel_location.values, bins=bin_edges
        )
        # norm_hist = this_cond_mod_hist
        # norm_hist = smooth(this_cond_mod_hist) / smooth(this_cond_mod_hist_all)
        norm_hist = this_cond_mod_hist / this_cond_mod_hist_all
        norm_hist = smooth(norm_hist)
        all_norm_hist.append(norm_hist)
        plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=cond_names[this_cond])
    plt.legend(loc='upper left', fontsize=10)
    plt.axvline(x=1, c='gray', alpha=0.5)
    plt.axvline(x=2, c='gray', alpha=0.5)
    plt.xticks([0.5, 1.5, 2.5], ['5\' UTR', 'CDS', '3\' UTR'])
    plt.ylabel('$N_{S\geq50}$ / $N_{covered}$', fontsize=12)
    plt.ylim([0, 0.1])
    plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
ymax = np.round((np.max(all_norm_hist) // 0.05 + 1) * 0.05, 2)
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    plt.ylim([0, ymax])
plt.suptitle(f'{ds}', fontsize=20)
plt.savefig(os.path.join(img_out, f"profile_{ds}_{'_'.join(conditions)}.png"), bbox_inches='tight')