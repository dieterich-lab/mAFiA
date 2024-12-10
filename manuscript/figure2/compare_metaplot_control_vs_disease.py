import os
import pandas as pd
import numpy as np
from scipy.stats import kstest
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
#######################################################################
import matplotlib.pyplot as plt


def smooth(y, box_pts=10):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']

ds = 'HFpEF'
conditions = ['ctrl_merged', 'HFpEF_merged']

# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']


ds_colors = {
    'HFpEF': 'r',
    'Diet': 'b',
    'TAC': 'g'
}

cond_names = {this_cond: this_cond.rstrip('_merged') for this_cond in conditions}
cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['gray', ds_colors[ds]])}

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresholds = ['0.0', '50.0']

metaplot_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

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

ymax = 0.08
plt.figure(figsize=(6*cm, 3*cm))
for mod_ind, this_mod in enumerate(mods):
    ax = plt.subplot(1, 2, mod_ind+1)
    all_norm_hist = []
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
        all_norm_hist.append(norm_hist)
        norm_hist = smooth(norm_hist)
        plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=cond_names[this_cond])
    ks_dist, pval = kstest(*all_norm_hist)
    # plt.legend(loc='upper left')
    plt.text(0.95, 0.05, rf'$\Delta$KS={ks_dist:.3f}' + f'\np-value={pval:.3E}', ha='right', va='bottom', transform=ax.transAxes)
    plt.axvline(x=1, c='gray', alpha=0.5)
    plt.axvline(x=2, c='gray', alpha=0.5)
    plt.xticks([0.5, 1.5, 2.5], ['5\' UTR', 'CDS', '3\' UTR'])
    # plt.xticks([0.5, 1.5, 2.5], [])
    # plt.ylabel('$N_{S\geq50}$ / $N_{covered}$')
    # plt.title(rf'${{{dict_mod_display[this_mod]}}}$')
# ymax = np.round((np.max(all_norm_hist) // 0.05 + 1) * 0.05, 2)
yticks = np.linspace(0, ymax, 3)
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    # plt.ylim([0, ymax])
    if mod_ind==0:
        plt.yticks(yticks, yticks)
    else:
        plt.yticks(yticks, [])
# plt.suptitle(f'{ds}')
plt.savefig(os.path.join(img_out, f"profile_{ds}_{'_'.join(conditions)}.{FMT}"), **fig_kwargs)