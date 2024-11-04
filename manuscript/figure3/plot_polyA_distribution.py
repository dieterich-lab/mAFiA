import os
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import pysam
from tqdm import tqdm
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200,
                  # transparent=True
                  )
#######################################################################
import matplotlib.pyplot as plt


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

########################################################################################################################
### define datasets ####################################################################################################
########################################################################################################################
# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']

# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

ds = 'Diet'
conditions = ['WT_CD_merged', 'WT_WD_merged']

cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['b', 'r'])}

polyA_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA'
polyA_files = {this_cond: os.path.join(polyA_dir, f'{ds}_{this_cond}_polyA.txt') for this_cond in conditions}

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'
os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### overall polyA distribution #########################################################################################
########################################################################################################################
df_polyA = {this_cond: pd.read_csv(polyA_files[this_cond], sep='\t', names=['read_id', 'polyA_len'])
            for this_cond in conditions}
dict_polyA = {this_cond: {k: v for k, v in df_polyA[this_cond][['read_id', 'polyA_len']].values}
              for this_cond in conditions}

xmax = 500
ymax = 0.25
xticks = np.linspace(0, xmax, 5)
yticks = np.arange(0, ymax+0.01, 0.1)

bin_width = 1
num_bins = np.int64(xmax / bin_width)
bin_edges = np.linspace(0, xmax, num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

### linear only ###
plt.figure(figsize=(4*cm, 4*cm))
# this_cond = conditions[0]
for this_cond in conditions:
    # plt.hist(dict_polyA[this_cond].values(), range=[0, 500], bins=50, label=this_cond, alpha=0.5, density=True)
    this_hist, _ = np.histogram(list(dict_polyA[this_cond].values()), bins=bin_edges)
    num_reads = np.sum(this_hist)
    norm_hist = this_hist / num_reads
    plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=this_cond.rstrip('_merged') + f' ({num_reads})')
plt.legend()
plt.xlabel('polyA length (bps)')
plt.ylabel('Norm. Frequency')
plt.xticks(xticks)
# plt.yticks(yticks)
plt.xlim([0, xmax])
# plt.ylim([0, ymax])
plt.savefig(os.path.join(img_out, f'hist_polyA_len_{ds}.{FMT}'), **fig_kwargs)

### linear and log ###
# plt.figure(figsize=(8*cm, 3*cm))
# plt.subplot(1, 2, 1)
# for cond_ind, this_cond in enumerate(conditions):
#     # plt.hist(dict_polyA[this_cond].values(), range=[0, 500], bins=50, label=this_cond, alpha=0.5, density=True)
#     this_hist, _ = np.histogram(list(dict_polyA[this_cond].values()), bins=bin_edges)
#     norm_hist = this_hist / np.sum(this_hist)
#     plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=this_cond)
# # plt.legend(fontsize=10)
# # plt.xlabel('polyA length (bps)', fontsize=12)
# # plt.ylabel('Density', fontsize=12)
# plt.subplot(1, 2, 2)
# for cond_ind, this_cond in enumerate(conditions):
#     # plt.hist(dict_polyA[this_cond].values(), range=[0, 500], bins=50, label=this_cond, alpha=0.5, density=True)
#     this_hist, _ = np.histogram(list(dict_polyA[this_cond].values()), bins=bin_edges)
#     norm_hist = this_hist / np.sum(this_hist)
#     plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=this_cond)
#     plt.yscale('log')
# plt.minorticks_off()
# plt.yticks([10**-5, 10**-3, 10**-1])
# # plt.legend(fontsize=10)
# # plt.xlabel('polyA length (bps)', fontsize=12)
# # plt.ylabel('Density (log)', fontsize=12)
# # plt.suptitle(ds, fontsize=15)
# plt.savefig(os.path.join(img_out, f'hist_polyA_len_{ds}_{conditions[0]}_vs_{conditions[1]}.{FMT}'), **fig_kwargs)
