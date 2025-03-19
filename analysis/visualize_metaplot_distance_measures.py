import os
import pandas as pd
import numpy as np
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 6
mpl.rcParams['legend.fontsize'] = 4
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
#######################################################################
mpl.use('TkAgg')
import matplotlib.pyplot as plt


def smooth(y, box_pts=5):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


ds = 'WT'

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
mod_colors = {
    'm6A': 'r',
    'psi': 'b'
}

thresholds = ['0.0', '50.0']

metaplot_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/metaPlotR'
img_out = '/home/adrian/img_out/manuscript_bioinformatics_application_note'
os.makedirs(img_out, exist_ok=True)

mod_thresh_dist_measure = {
    this_mod: {} for this_mod in mods
}

for this_mod in mods:
    for this_thresh in thresholds:
        this_dist_measure_file = os.path.join(metaplot_dir, f"{ds}_{this_mod}_modRatio{this_thresh}.dist.measures.txt")
        mod_thresh_dist_measure[this_mod][this_thresh] = pd.read_csv(this_dist_measure_file, sep='\t')

hist_range = [0, 3]
hist_num_bins = 60
bin_edges = np.linspace(*hist_range, hist_num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

plt.figure(figsize=(8*cm, 3*cm))
all_norm_hist = []
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    this_mod_hist, _ = np.histogram(
        mod_thresh_dist_measure[this_mod][thresholds[1]].rel_location.values, bins=bin_edges
    )
    this_mod_hist_all, _ = np.histogram(
        mod_thresh_dist_measure[this_mod]['0.0'].rel_location.values, bins=bin_edges
    )
    # norm_hist = this_cond_mod_hist
    # norm_hist = smooth(this_cond_mod_hist) / smooth(this_cond_mod_hist_all)
    norm_hist = this_mod_hist / this_mod_hist_all
    norm_hist = smooth(norm_hist)
    all_norm_hist.append(norm_hist)
    plt.plot(bin_centers, norm_hist, c=mod_colors[this_mod])
    # plt.legend(loc='upper left', fontsize=10)
    plt.axvline(x=1, c='gray', alpha=0.5)
    plt.axvline(x=2, c='gray', alpha=0.5)
    plt.xticks([0.5, 1.5, 2.5], ['5\' UTR', 'CDS', '3\' UTR'])
    # plt.ylabel('$N_{S\geq50}$ / $N_{covered}$', fontsize=12)
    # plt.ylim([0, 0.1])
    # plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
ymax = np.round((np.max(all_norm_hist) // 0.05 + 1) * 0.05, 2)
yticks = np.linspace(0, ymax, 4)
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    plt.ylim([0, ymax])
    plt.yticks(yticks)
# plt.suptitle(f'{ds}', fontsize=20)
plt.savefig(os.path.join(img_out, f"metagene_profile_{ds}.{FMT}"), **fig_kwargs)
