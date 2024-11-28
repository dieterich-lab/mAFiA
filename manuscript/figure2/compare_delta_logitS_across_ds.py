import os
import pandas as pd
from functools import reduce
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
import matplotlib.tri as tri
import scipy

in_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_logit'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

# ds = ['TAC', 'HFpEF', 'Diet']
ds = ['TAC', 'HFpEF']

ds_display = {
    'TAC': 'TAC /\nSHAM',
    'HFpEF': 'HFpEF /\nCtrl',
    # 'Diet': 'WD /\nCD'
}
ds_color = {
    'TAC': 'g',
    'HFpEF': 'r',
    'Diet': 'b'
}
mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresh_pval = 0.05
num_gene_mods = 50

dfs = {}
for this_ds in ds:
    this_df = pd.read_csv(os.path.join(in_dir, f'delta_logitS_{this_ds}.tsv'), sep='\t')
    this_df = this_df[this_df.pval < thresh_pval]
    this_df.rename(columns={
        'num_sites': f'num_sites_{this_ds}',
        'delta_logit': f'delta_logit_{this_ds}',
        'pval': f'pval_{this_ds}'
    }, inplace=True)
    dfs[this_ds] = this_df

df_merged = reduce(lambda left, right: pd.merge(left, right, on=['gene', 'mod']), list(dfs.values()))

### split by mod and ds ###
vmax = 0.5
# fig = plt.figure()
fig, axes = plt.subplots(figsize=(10*cm, 12*cm), nrows=1, ncols=4)
for mod_ind, this_mod in enumerate(mods):
    # mod_df_merged = df_merged[df_merged['mod'] == this_mod]
    for ds_ind, this_ds in enumerate(ds):
        df_ds = dfs[this_ds]
        mod_df_merged = df_ds[df_ds['mod'] == this_mod]
        ds_mod_df_merged = mod_df_merged.sort_values(f'num_sites_{this_ds}', ascending=False)[:num_gene_mods]
        # ds_mod_df_merged = mod_df_merged.sort_values(f'pval_{this_ds}', ascending=True)[:num_gene_mods]
        ds_mod_df_merged.sort_values(f'delta_logit_{this_ds}', ascending=False, inplace=True)

        vec_gene = ds_mod_df_merged['gene'].values
        vec_delta_logit = ds_mod_df_merged[f'delta_logit_{this_ds}'].values

        subplot_ind = mod_ind*2+ds_ind+1
        plt.subplot(1, 4, subplot_ind)
        im = plt.imshow(vec_delta_logit[:, np.newaxis], aspect=0.5, origin='upper', vmin=-vmax, vmax=vmax, cmap='seismic')
        plt.xticks([0], [this_ds])
        plt.yticks(np.arange(num_gene_mods), vec_gene)
        plt.gca().xaxis.tick_top()

fig.subplots_adjust(top=1, bottom=0.05)
cbar_ax = fig.add_axes([0.25, 0.01, 0.5, 0.02])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_ticks(np.linspace(-vmax, vmax, 5))
plt.savefig(os.path.join(img_out, f'heatmap_delta_logit_S_split.{FMT}'), **fig_kwargs)

########################################################################################################################

vec_gene_mod = [f'{gene}, ${dict_mod_display[mod]}$' for gene, mod in df_merged[['gene', 'mod']].values]
mat_delta_logit = df_merged[[f'delta_logit_{this_ds}' for this_ds in ds]].values

### select rows ###
# vec_variance = np.sum(np.abs(mat_delta_logit), axis=1)
# vec_variance = np.sum(mat_delta_logit**2, axis=1)
vec_variance = df_merged[[f'num_sites_{this_ds}' for this_ds in ds]].values.sum(axis=1)

sel_indices = np.argpartition(vec_variance, -num_gene_mods)[-num_gene_mods:]

for row in [vec_gene_mod[ind] for ind in sel_indices]:
    print(row)

sel_gene_mod = [vec_gene_mod[this_ind] for this_ind in sel_indices]
sel_mat_delta_logit = mat_delta_logit[sel_indices]

# sorted_indices = np.lexsort((
#         sel_mat_delta_logit[:, 0],
#         # np.round(sel_mat_delta_logit[:, 1]),
#         # np.round(sel_mat_delta_logit[:, 2])
#     ))
sorted_indices = np.argsort(sel_gene_mod)[::-1]
sorted_mat_delta_logit = sel_mat_delta_logit[sorted_indices]
sorted_gene_mod = np.array(sel_gene_mod)[sorted_indices]

vmax = 0.5
plt.figure(figsize=(4*cm, 15*cm))
plt.imshow(sorted_mat_delta_logit, aspect='auto', origin='lower', vmin=-vmax, vmax=vmax, cmap='seismic')
plt.xticks(np.arange(len(ds)), list(ds_display.values()))
plt.yticks(np.arange(num_gene_mods), sorted_gene_mod)
plt.gca().xaxis.tick_top()
cbar = plt.colorbar()
cbar.set_ticks(np.linspace(-vmax, vmax, 5))
plt.savefig(os.path.join(img_out, f'heatmap_delta_logit_S_all_ds.{FMT}'), **fig_kwargs)

### scatter plots in (m6A, psi) space ###
ds_vec_mods = {this_ds: [] for this_ds in ds}
for this_ds in ds:
    this_df = dfs[this_ds]
    unique_genes = this_df['gene'].unique()
    for this_gene in unique_genes:
        sub_df = this_df[
            (this_df['gene'] == this_gene)
            # * (this_df[f'pval_{this_ds}'] < 0.001)
            ]
        if (len(sub_df['mod']) == 2) and set(mods).issubset(set(sub_df['mod'].values)):
            sub_df.sort_values('mod', inplace=True)
            ds_vec_mods[this_ds].append(sub_df[f'delta_logit_{this_ds}'].values)

ds_cmap = {
    'TAC': 'Greens',
    'HFpEF': 'Reds',
    'Diet': 'Blues'
}

xmax = 0.7
zoom_factor = 2
num_bins = 20
levels = np.linspace(0, 2, 5)
ticks = np.linspace(-xmax, xmax, 3)
plt.figure(figsize=(15*cm, 4*cm))
for ds_ind, this_ds in enumerate(ds):
    plt.subplot(1, 3, ds_ind+1)
    plt.axvline(x=0, c='gray', ls='--', alpha=0.5)
    plt.axhline(y=0, c='gray', ls='--', alpha=0.5)
    vec_x, vec_y = np.vstack(ds_vec_mods[this_ds]).T
    # plt.scatter(vec_x, vec_y, s=1, c=ds_color[this_ds], label=ds_display[this_ds])
    mat_z, bin_y, bin_x = np.histogram2d(vec_x, vec_y, range=[[-xmax, xmax], [-xmax, xmax]], bins=num_bins)
    mat_z = scipy.ndimage.zoom(mat_z, zoom_factor)
    mat_z[mat_z < 10] = 0
    delta_x = bin_x[1] - bin_x[0]
    mat_z = mat_z / (np.sum(mat_z) * delta_x**2)
    # center_x = 0.5 * (bin_x[1:] + bin_x[:-1])
    zoom_bin_x = np.linspace(-xmax, xmax, num_bins*zoom_factor+1)
    zoom_center_x = 0.5 * (zoom_bin_x[1:] + zoom_bin_x[:-1])
    mat_x, mat_y = np.meshgrid(zoom_center_x, zoom_center_x)
    plt.contourf(mat_x, mat_y, mat_z, levels=levels, cmap=ds_cmap[this_ds], vmin=0, vmax=2)
    # plt.contourf(mat_z, levels=3, cmap=ds_cmap[this_ds], vmin=20, alpha=0.5)
    plt.xticks(ticks)
    plt.yticks(ticks)
    # if ds_ind == 0:
    #     plt.yticks(ticks)
    # else:
    #     plt.yticks(ticks, [])
    plt.xlim([-xmax, xmax])
    plt.ylim([-xmax, xmax])
    plt.xlabel('delta logit S, m6A')
    plt.ylabel('delta logit S, psi')
    # plt.colorbar(orientation='horizontal', location='top', ticks=[0, 1, 2])
    plt.colorbar(ticks=[0, 1, 2])
plt.savefig(os.path.join(img_out, f'contour_delta_logit_S_all_ds.{FMT}'), **fig_kwargs)
