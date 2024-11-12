import os
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
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


def get_logit(vec_modRatio):
    return np.log2(vec_modRatio / (100 - vec_modRatio))


bed7_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

thresh_conf = 80
thresh_cov = 50

ds = ['TAC', 'HFpEF']
conditions = ['TAC_merged', 'HFpEF_merged']

dfs = []
for this_ds, this_cond in zip(ds, conditions):
    res_file = os.path.join(res_dir, f'{this_ds}/{this_cond}/chrALL.mAFiA.sites.bed')
    df_res = pd.read_csv(res_file, sep='\t', dtype={'chrom': str})
    df_res = df_res[
        (df_res['confidence'] >= thresh_conf)
        * (df_res['coverage'] >= thresh_cov)
    ]
    df_res.rename(columns={
        'coverage': f'coverage_{this_ds}',
        'confidence': f'confiednece_{this_ds}',
        'modRatio': f'modRatio_{this_ds}'
    }, inplace=True)
    dfs.append(df_res)

df_merged = pd.merge(*dfs, on=bed7_fields)

ticks = np.linspace(0, 100, 3)
num_bins = 20
# plt.figure(figsize=(8*cm, 8*cm))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(2, 2, mod_ind+1)
#     mod_df = df_merged[df_merged['name'] == this_mod]
#     vec_x, vec_y = mod_df[[f'modRatio_{this_ds}' for this_ds in ds]].values.T
#
#     delta_logit_S = get_logit(vec_y) - get_logit(vec_x)
#
#     plt.scatter(vec_x, vec_y, s=1, rasterized=True)
#
#     # grid_z, bin_x, bin_y = np.histogram2d(vec_x, vec_y, bins=num_bins, range=[[0, 100], [0, 100]])
#     # grid_z = grid_z / grid_z.sum()
#     # log_grid_z = np.log10(grid_z+1)
#     # bin_centers = 0.5 * (bin_x[1:] + bin_x[:-1])
#     # grid_x, gird_y = np.meshgrid(bin_centers, bin_centers)
#     # plt.contourf(grid_x, gird_y, grid_z, levels=10, vmin=0, vmax=0.1)
#
#     plt.plot([0, 100], [0, 100], c='r', ls='--')
#     plt.title(rf'{len(vec_x)} common ${{{dict_mod_display[this_mod]}}}$ sites')
#     plt.xticks(ticks)
#     plt.yticks(ticks)
#     plt.xlim([-1, 101])
#     plt.ylim([-1, 101])
#
#     bin_max = 3
#     xticks = np.linspace(-bin_max, bin_max, 3)
#     plt.subplot(2, 2, mod_ind+3)
#     plt.hist(delta_logit_S, bins=60, range=[-bin_max, bin_max])
#     plt.axvline(x=0, c='r', ls='--')
#     plt.xticks(xticks)
#     plt.xlim([-bin_max, bin_max])
#
# plt.savefig(os.path.join(img_out, f"delta_logit_S_{'_'.join(ds)}.{FMT}"), **fig_kwargs)

bin_max = 3
xticks = np.linspace(-bin_max, bin_max, 3)
plt.figure(figsize=(8*cm, 3*cm))
for mod_ind, this_mod in enumerate(mods):
    mod_df = df_merged[df_merged['name'] == this_mod]
    vec_x, vec_y = mod_df[[f'modRatio_{this_ds}' for this_ds in ds]].values.T
    delta_logit_S = get_logit(vec_y) - get_logit(vec_x)

    plt.subplot(1, 2, mod_ind+1)
    plt.hist(delta_logit_S, bins=60, range=[-bin_max, bin_max])
    plt.axvline(x=0, c='r', ls='--')
    plt.xticks(xticks)
    plt.xlim([-bin_max, bin_max])
    plt.title(f'N = {len(vec_x)}')

plt.savefig(os.path.join(img_out, f"hist_delta_logit_S_{'_'.join(ds)}.{FMT}"), **fig_kwargs)


