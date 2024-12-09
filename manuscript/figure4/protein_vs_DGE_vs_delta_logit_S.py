import os
from glob import glob
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import numpy as np
# from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
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
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200, transparent=True)
#######################################################################


def plot_distribution(vec_mod, vec_mask, thresh_mask, mod_name, bin_max, ymax=None):
    # bin_max = 0.5
    # bin_width = 0.3
    # num_bins = int(2 * bin_max / bin_width)
    num_bins = 5
    bin_edges = np.round(np.linspace(-bin_max, bin_max, num_bins + 1), 2)
    bin_width = np.round(bin_edges[1] - bin_edges[0], 2)
    bin_centers = np.round(0.5 * (bin_edges[1:] + bin_edges[:-1]), 2)

    hist_pos, _ = np.histogram(vec_mod[vec_mask > thresh_mask], range=[-bin_max, bin_max], bins=num_bins)
    num_pos = np.sum(hist_pos)
    norm_pos = hist_pos / num_pos
    hist_neg, _ = np.histogram(vec_mod[vec_mask < -thresh_mask], range=[-bin_max, bin_max], bins=num_bins)
    num_neg = np.sum(hist_neg)
    norm_neg = hist_neg / num_neg
    plt.bar(bin_centers+bin_width*0.2, norm_pos, label=f'Up-regulated ({num_pos})', color='r', width=bin_width*0.4)
    plt.bar(bin_centers-bin_width*0.2, norm_neg, label=f'Down-regulated ({num_neg})', color='b', width=bin_width*0.4)
    # plt.plot(bin_centers, norm_pos, label=f'Up-regulated ({num_pos})', color='r')
    # plt.plot(bin_centers, norm_neg, label=f'Down-regulated ({num_neg})', color='b')
    plt.legend(loc='upper left')
    # plt.xlabel(f'Delta logit S, {mod_name}')
    # plt.ylabel('Prob(DGE)')
    # plt.ylim([0, ymax])
    plt.xticks(bin_centers, bin_centers)
    if ymax is not None:
        plt.yticks(np.round(np.linspace(0, ymax, 3), 2))
    return norm_pos, norm_neg, bin_centers, bin_width


def plot_prob_diff(norm_pos, norm_neg, mod_name, bin_centers, bin_width, ymax=None):
    hist_diff = (norm_pos - norm_neg) / (norm_pos + norm_neg + 0.000001)
    hist_pos = hist_diff * (hist_diff>=0)
    hist_neg = hist_diff * (hist_diff<0)
    plt.bar(bin_centers, hist_pos, color='r', width=bin_width*0.8)
    plt.bar(bin_centers, hist_neg, color='b', width=bin_width*0.8)
    plt.axhline(y=0, c='gray', ls='--')
    # plt.xlabel(f'Delta logit S, {mod_name}')
    # plt.ylabel('Prob(DGE | mod)')
    plt.xticks(bin_centers, bin_centers)
    if ymax is not None:
        plt.ylim([-ymax, ymax])
        plt.yticks(np.round(np.linspace(-ymax, ymax, 5), 2))


img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
abundance_dir = os.path.join(res_dir, 'DTU')
logit_dir = os.path.join(res_dir, 'transcript_logit')
gene_count_dir = os.path.join(res_dir, 'gene_counts')
protein_dir = os.path.join(res_dir, 'protein_abundance')
polyA_dir = os.path.join(res_dir, 'polyA')

mods = ['m6A', 'psi']

thresh_fdr = 0.1
thresh_pval = 0.01

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']
# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

df_abundance = pd.read_excel(
    glob(os.path.join(abundance_dir, f'DTU_{ds}_*.xlsx'))[0],
    'DGE'
)
df_abundance.rename(columns={'mgi_symbol': 'gene'}, inplace=True)

# df_abundance = df_abundance[df_abundance['FDR'] < thresh_fdr]
df_abundance = df_abundance[
    ['BambuGene' not in this_id for this_id in df_abundance['ID']]
]

### add gene counts ###
dfs_gene_count = []
for this_cond in conditions:
    dfs_gene_count.append(
        pd.read_csv(os.path.join(gene_count_dir, f'{ds}_{this_cond}.counts'), sep='\t', names=['ID', 'gene_count'])
    )
df_gene_count = pd.merge(*dfs_gene_count, on=['ID'], suffixes=[f'_{this_cond}' for this_cond in conditions])

df_append = df_gene_count.set_index('ID').loc[df_abundance['ID']][[
    f'gene_count_{this_cond}' for this_cond in conditions
]].reset_index()
df_abundance = pd.merge(df_abundance, df_append, on=['ID'])
# out_abundance_filename = os.path.join(abundance_dir, f'DGE_{ds}_with_gene_count.tsv')
# df_abundance.to_csv(out_abundance_filename, sep='\t', index=False)
df_abundance.sort_values(f'gene_count_{conditions[0]}', ascending=False, inplace=True)
df_abundance.rename(columns={'logFC': 'log2fc_gene'}, inplace=True)


### volcano plot ###
xmax = 4
ymax = 1
num_genes_display = 10
margin = 0.95
xticks = np.linspace(-xmax, xmax, 5)
yticks = np.linspace(0, 1, 5)
plt.figure(figsize=(4*cm, 4*cm))
vec_log2fc_gene = df_abundance['log2fc_gene'].values
vec_inv_fdr = 1.0 - df_abundance['FDR'].values
mask = vec_inv_fdr > (1.0 - thresh_fdr)
mask_pos = mask * (vec_log2fc_gene > 0)
mask_neg = mask * (vec_log2fc_gene < 0)
gene_pos = df_abundance['gene'].values[mask_pos]
gene_neg = df_abundance['gene'].values[mask_neg]
plt.scatter(vec_log2fc_gene[~mask], vec_inv_fdr[~mask], s=0.5, c='gray', alpha=0.5, edgecolors='none', rasterized=True)
plt.scatter(vec_log2fc_gene[mask_pos], vec_inv_fdr[mask_pos], s=0.5, c='red', edgecolors='none', rasterized=True)
plt.scatter(vec_log2fc_gene[mask_neg], vec_inv_fdr[mask_neg], s=0.5, c='blue', edgecolors='none', rasterized=True)
plt.text(xmax*margin, ymax*margin, '\n'.join(gene_pos[:num_genes_display]), ha='right', va='top', c='red')
plt.text(-xmax*margin, ymax*margin, '\n'.join(gene_neg[:num_genes_display]), ha='left', va='top', c='blue')
plt.text(xmax*margin, ymax*margin+0.1, f'{len(gene_pos)} up', ha='right', va='bottom', c='red')
plt.text(-xmax*margin, ymax*margin+0.1, f'{len(gene_neg)} down', ha='left', va='bottom', c='blue')
# plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]) + f'\n...\nTotal: {len(gene_pos)}', ha='right', va='top', c='red')
# plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]) + f'\n...\nTotal: {len(gene_neg)}', ha='left', va='top', c='blue')
# plt.title(f'{len(gene_pos)} up, {len(gene_neg)} down')
plt.xlim([-xmax, xmax])
plt.ylim([0, ymax + 0.01])
plt.xticks(xticks)
plt.yticks(yticks)
# plt.xlabel('log2fc gene expression')
# plt.ylabel('1 - FDR')
plt.savefig(os.path.join(img_out, f'volcano_DGE_{ds}_FDR{thresh_fdr}.{FMT}'), **fig_kwargs)


### protein abundance ###
df_protein = pd.read_csv(
    os.path.join(protein_dir, 'TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'),
    sep='\t'
)
df_protein = df_protein[df_protein['tp'] == 'Combined']
# df_protein.sort_values('p.adj', inplace=True)
df_protein.rename(columns={'SYMBOL': 'gene'}, inplace=True)

log2fc_protein = []
for this_gene in df_abundance['gene']:
    sub_df_protein = df_protein[df_protein['gene'] == this_gene]
    if len(sub_df_protein) == 1:
        log2fc_protein.append(sub_df_protein['logFC'].values[0])
    else:
        log2fc_protein.append(np.nan)
df_abundance['log2fc_protein'] = log2fc_protein

df_valid = df_abundance[
    ~np.isnan(df_abundance['log2fc_gene'])
    * ~np.isnan(df_abundance['log2fc_protein'])
]
mask_up = ((df_valid['FDR'].values < thresh_fdr) * (df_valid['log2fc_gene'] > 0)).values
mask_down = ((df_valid['FDR'].values < thresh_fdr) * (df_valid['log2fc_gene'] < 0)).values


### scatter plot protein vs gene ###
def fit_func(m, b, x):
    return b+m*x

xmax = 3.2
ymax = 1.1
# vec_x = df_valid['log2fc_gene'].values
# vec_y = df_valid['log2fc_protein'].values
# mask_fit = (np.abs(vec_x) < xmax) * (np.abs(vec_y) < ymax)
# popt, pcov = curve_fit(fit_func, vec_x[mask_fit], vec_y[mask_fit])

xticks = np.linspace(-int(xmax), int(xmax), 5)
yticks = np.linspace(-int(ymax), int(ymax), 5)
plt.figure(figsize=(4*cm, 4*cm))
plt.scatter(df_valid['log2fc_gene'], df_valid['log2fc_protein'],
            s=0.5, edgecolor='none', c='gray', alpha=0.5, rasterized=True)
plt.scatter(df_valid['log2fc_gene'][mask_up], df_valid['log2fc_protein'][mask_up],
            s=1, c='red', rasterized=True)
plt.scatter(df_valid['log2fc_gene'][mask_down], df_valid['log2fc_protein'][mask_down],
            s=1, c='blue', rasterized=True)
# plt.plot(xticks, fit_func(xticks, *popt), ls='--', alpha=0.5)
plt.xlim([-xmax, xmax])
plt.ylim([-ymax, ymax])
plt.xticks(xticks)
plt.yticks(yticks)

plt.savefig(os.path.join(img_out, f'scatter_DGE_vs_protein_{ds}_FDR{thresh_fdr}.{FMT}'), **fig_kwargs)


### add delta logit S ###
df_logit = pd.read_csv(
    os.path.join(logit_dir, f'delta_logitS_{ds}.tsv'),
    sep='\t'
)

gene_mod_logit = []
for this_gene in df_abundance['gene']:
    sub_df_logit = df_logit[df_logit['gene'] == this_gene]
    this_gene_mod_logit = [this_gene]
    if len(sub_df_logit):
        for this_mod in mods:
            mod_logit = sub_df_logit[sub_df_logit['mod'] == this_mod]
            if len(mod_logit) == 1:
                this_gene_mod_logit.append(mod_logit['delta_logit'].values[0])
            else:
                this_gene_mod_logit.append(np.nan)
    else:
        this_gene_mod_logit = [this_gene, np.nan, np.nan]
    gene_mod_logit.append(this_gene_mod_logit)
df_mod_logit = pd.DataFrame(gene_mod_logit, columns=['gene'] + [f'delta_logit_S_{this_mod}' for this_mod in mods])
df_abundance_mod = pd.merge(df_abundance, df_mod_logit, on=['gene'])


### (m6A, psi) scatter ###
from scipy import interpolate
from scipy.ndimage import zoom

xymax = 0.5

def get_binned_z(in_vec_x, in_vec_y, in_vec_z, bin_max=0.5, num_bins=10, method='sum'):
    bin_edges = np.linspace(-bin_max, bin_max, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    grid_z = np.zeros((num_bins, num_bins))
    for bin_i in range(num_bins):
        i_start = bin_edges[bin_i]
        i_end = bin_edges[bin_i+1]
        i_mask = (in_vec_y >= i_start) * (in_vec_y < i_end)
        for bin_j in range(num_bins):
            j_start = bin_edges[bin_j]
            j_end = bin_edges[bin_j+1]
            j_mask = (in_vec_x >= j_start) * (in_vec_x < j_end)
            binned_z = in_vec_z[i_mask*j_mask]
            if len(binned_z):
                if method == 'sum':
                    grid_z[bin_i, bin_j] = np.sum(binned_z)
                elif method == 'mean':
                    grid_z[bin_i, bin_j] = np.nanmean(binned_z)
    grid_x, grid_y = np.meshgrid(bin_centers, bin_centers)
    return grid_x, grid_y, grid_z


df_up = df_abundance_mod[
    (df_abundance_mod['FDR'] < thresh_fdr)
    * (df_abundance_mod['log2fc_gene'] > 0)
]
vec_up_x, vec_up_y, vec_up_z = df_up[['delta_logit_S_m6A', 'delta_logit_S_psi', 'log2fc_gene']].values.T
up_grid_x, up_grid_y, up_grid_z = get_binned_z(vec_up_x, vec_up_y, vec_up_z, xymax)

df_down = df_abundance_mod[
    (df_abundance_mod['FDR'] < thresh_fdr)
    * (df_abundance_mod['log2fc_gene'] < 0)
]
vec_down_x, vec_down_y, vec_down_z = df_down[['delta_logit_S_m6A', 'delta_logit_S_psi', 'log2fc_gene']].values.T
down_grid_x, down_grid_y, down_grid_z = get_binned_z(vec_down_x, vec_down_y, vec_down_z, xymax)

up_grid_z = up_grid_z / up_grid_z.sum()
down_grid_z = down_grid_z / down_grid_z.sum()

zoom_factor = 3
up_grid_x_zoom = zoom(up_grid_x, zoom_factor)
up_grid_y_zoom = zoom(up_grid_y, zoom_factor)
up_grid_z_zoom = zoom(up_grid_z, zoom_factor)
down_grid_x_zoom = zoom(down_grid_x, zoom_factor)
down_grid_y_zoom = zoom(down_grid_y, zoom_factor)
down_grid_z_zoom = zoom(down_grid_z, zoom_factor)

up_grid_z_zoom[up_grid_z_zoom < np.quantile(up_grid_z_zoom, 0.25)] = 0
down_grid_z_zoom[down_grid_z_zoom < np.quantile(down_grid_z_zoom, 0.25)] = 0

vmax = 0.1
n_levels = 3
ticks = np.linspace(-xymax, xymax, 3)

plt.figure(figsize=(10*cm, 4*cm))
plt.subplot(1, 2, 1)
plt.contourf(up_grid_x_zoom, up_grid_y_zoom, up_grid_z_zoom,
             vmin=0, vmax=vmax, levels=n_levels, cmap='Reds', alpha=0.5)
plt.axvline(x=0, c='gray', ls='--')
plt.axhline(y=0, c='gray', ls='--')
plt.xticks(ticks)
plt.yticks(ticks)
cbar = plt.colorbar()
cbar.set_ticks(np.linspace(0, vmax, 3))
plt.subplot(1, 2, 2)
plt.contourf(down_grid_x_zoom, down_grid_y_zoom, down_grid_z_zoom,
             vmin=0, vmax=vmax, levels=n_levels, cmap='Blues', alpha=0.5)
plt.axvline(x=0, c='gray', ls='--')
plt.axhline(y=0, c='gray', ls='--')
plt.xticks(ticks)
plt.yticks(ticks)
cbar = plt.colorbar()
cbar.set_ticks(np.linspace(0, vmax, 3))

plt.savefig(os.path.join(img_out, f'contour_DGE_in_mod_plane_{ds}_FDR{thresh_fdr}.{FMT}'), **fig_kwargs)

### polyA ###
# df_polyA = pd.read_csv(
#     os.path.join(polyA_dir, f'gene_polyA_log2fc_pval_{ds}.tsv'),
#     sep='\t'
# )
#
# gene_polyA = []
# for this_gene in df_abundance['gene']:
#     sub_df_polyA = df_polyA[df_polyA['gene'] == this_gene]
#     if len(sub_df_polyA) == 1:
#         gene_polyA.append([this_gene] + list(sub_df_polyA[['log2fc', 'pval']].values[0]))
#     else:
#         gene_polyA.append([this_gene, np.nan, np.nan])
# df_gene_polyA = pd.DataFrame(gene_polyA, columns=['gene', 'log2fc_polyA', 'pval_polyA'])
# df_abundance_mod_polyA = pd.merge(df_abundance_mod, df_gene_polyA, on=['gene'])
#
# df_polyA_up = df_abundance_mod_polyA[
#     (df_abundance_mod_polyA['log2fc_polyA'] > 0.01)
#     * (df_abundance_mod_polyA['pval_polyA'] < thresh_pval)
# ]
# vec_polyA_up_x, vec_polyA_up_y, vec_polyA_up_z = df_polyA_up[['delta_logit_S_m6A', 'delta_logit_S_psi', 'log2fc_polyA']].values.T
# polyA_up_grid_x, polyA_up_grid_y, polyA_up_grid_z = get_binned_z(vec_polyA_up_x, vec_polyA_up_y, vec_polyA_up_z,
#                                                                  xymax, method='mean')
#
# zoom_factor_polyA = 3
# polyA_up_grid_x_zoom = zoom(polyA_up_grid_x, zoom_factor_polyA)
# polyA_up_grid_y_zoom = zoom(polyA_up_grid_y, zoom_factor_polyA)
# polyA_up_grid_z_zoom = zoom(polyA_up_grid_z, zoom_factor_polyA)
#
# plt.contourf(polyA_up_grid_x_zoom, polyA_up_grid_y_zoom, polyA_up_grid_z_zoom,
#              vmin=0.1, vmax=0.8, levels=n_levels, cmap='Reds', alpha=0.5)
#
#
# df_polyA_down = df_abundance_mod_polyA[
#     (df_abundance_mod_polyA['log2fc_polyA'] < -0.01)
#     * (df_abundance_mod_polyA['pval_polyA'] < thresh_pval)
# ]
# vec_polyA_down_x, vec_polyA_down_y, vec_polyA_down_z = df_polyA_down[['delta_logit_S_m6A', 'delta_logit_S_psi', 'log2fc_polyA']].values.T
# polyA_down_grid_x, polyA_down_grid_y, polyA_down_grid_z = get_binned_z(vec_polyA_down_x, vec_polyA_down_y, vec_polyA_down_z,
#                                                                  xymax, method='mean')
#
# polyA_down_grid_x_zoom = zoom(polyA_down_grid_x, zoom_factor_polyA)
# polyA_down_grid_y_zoom = zoom(polyA_down_grid_y, zoom_factor_polyA)
# polyA_down_grid_z_zoom = zoom(polyA_down_grid_z, zoom_factor_polyA)
#
# plt.contourf(polyA_down_grid_x_zoom, polyA_down_grid_y_zoom, polyA_down_grid_z_zoom,
#              vmin=-0.8, vmax=0, levels=n_levels, cmap='Blues', alpha=0.5)


# plt.pcolormesh(up_grid_x, up_grid_y, up_grid_z, vmin=-vmax, vmax=vmax, cmap='bwr')
# plt.pcolormesh(down_grid_x, down_grid_y, down_grid_z, vmin=-vmax, vmax=vmax, cmap='bwr')



# thresh_coverage = 20
# num_sufficient_coverage = (df_abundance[
#                                [f'gene_count_{this_cond}' for this_cond in conditions]
#                            ].values >= thresh_coverage).all(axis=1).sum()
# with open(out_abundance_filename, 'a') as fout:
#     fout.write(f'\n\n#{ds}: {num_sufficient_coverage} / {len(df_abundance)} with coverage >= {thresh_coverage}')
#

### volcano plot ###
# xmax = 4
# ymax = 20
# xticks = np.linspace(-xmax, xmax, 5)
# yticks = np.linspace(0, ymax, 5)
#
# plt.figure(figsize=(8*cm, 4*cm))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     df_mod = df_logit[df_logit['mod'] == this_mod]
#     df_mod.sort_values('pval', inplace=True)
#     vec_delta_logit = df_mod['delta_logit'].values
#     vec_neg_log_pval = -np.log10(df_mod['pval'].values)
#     mask = vec_neg_log_pval > -np.log10(thresh_pval)
#     mask_pos = mask * (vec_delta_logit > 0)
#     mask_neg = mask * (vec_delta_logit < 0)
#     gene_pos = df_mod['gene'].values[mask_pos]
#     gene_neg = df_mod['gene'].values[mask_neg]
#     plt.scatter(vec_delta_logit[~mask], vec_neg_log_pval[~mask], s=0.5, c='gray', alpha=0.5, edgecolors='none', rasterized=True)
#     plt.scatter(vec_delta_logit[mask_pos], vec_neg_log_pval[mask_pos], s=0.5, c='red', edgecolors='none', rasterized=True)
#     plt.scatter(vec_delta_logit[mask_neg], vec_neg_log_pval[mask_neg], s=0.5, c='blue', edgecolors='none', rasterized=True)
#     plt.text(xmax*0.9, ymax*0.9, '\n'.join(gene_pos[:10]), ha='right', va='top', c='red')
#     plt.text(-xmax*0.9, ymax*0.9, '\n'.join(gene_neg[:10]), ha='left', va='top', c='blue')
#     # plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]) + f'\n...\nTotal: {len(gene_pos)}', ha='right', va='top', c='red')
#     # plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]) + f'\n...\nTotal: {len(gene_neg)}', ha='left', va='top', c='blue')
#     # plt.title(f'{len(gene_pos)} up, {len(gene_neg)} down')
#     plt.xlim([-xmax, xmax])
#     plt.ylim([0, ymax])
#     plt.xticks(xticks)
#     plt.yticks(yticks)
# # plt.xlabel('log2fc gene expression')
# # plt.ylabel('1 - FDR')
# plt.savefig(os.path.join(img_out, f'volcano_delta_logit_{ds}_pval{thresh_pval}.{FMT}'), **fig_kwargs)
#
# # df_logit = df_logit[df_logit['pval'] < thresh_pval]
#
# gene_dge_m6A_psi = []
# for this_gene, this_dge in df_abundance[['gene', 'logFC']].values:
#     sub_df_logit = df_logit[df_logit['gene'] == this_gene]
#     if 'm6A' in sub_df_logit['mod'].values:
#         this_m6A = sub_df_logit[sub_df_logit['mod'] == 'm6A']['delta_logit'].values[0]
#     else:
#         this_m6A = np.nan
#     if 'psi' in sub_df_logit['mod'].values:
#         this_psi = sub_df_logit[sub_df_logit['mod'] == 'psi']['delta_logit'].values[0]
#     else:
#         this_psi = np.nan
#     gene_dge_m6A_psi.append(
#         (
#             this_gene,
#             this_dge,
#             this_m6A,
#             this_psi
#         )
#     )

# vec_gene, vec_dge, vec_m6A, vec_psi = np.vstack(gene_dge_m6A_psi).T
# vec_dge = np.float64(vec_dge)
# vec_m6A = np.float64(vec_m6A)
# vec_psi = np.float64(vec_psi)
#
# thresh_dge = 0
#
# plt.figure(figsize=(8*cm, 8*cm))
# plt.subplot(2, 2, 1)
# norm_pos_m6A, norm_neg_m6A, centers, width = plot_distribution(vec_m6A, vec_dge, thresh_dge, 'm6A',
#                                                                bin_max=0.5, ymax=0.7)
# plt.subplot(2, 2, 2)
# plot_prob_diff(norm_pos_m6A, norm_neg_m6A, 'm6A', centers, width, ymax=0.6)
#
# plt.subplot(2, 2, 3)
# norm_pos_psi, norm_neg_psi, centers, width = plot_distribution(vec_psi, vec_dge, thresh_dge, 'psi',
#                                                                bin_max=0.5, ymax=0.7)
# plt.subplot(2, 2, 4)
# plot_prob_diff(norm_pos_psi, norm_neg_psi, 'psi', centers, width, ymax=0.6)
#
# # plt.suptitle(ds)
# plt.savefig(os.path.join(img_out, f'hist_DGE_vs_delta_logit_S_{ds}_FDR{thresh_fdr}_pval{thresh_pval}.{FMT}'),
#             **fig_kwargs)

### volcano plot ###
# xmax = 2
# ymax = 8
# xticks = np.linspace(-xmax, xmax, 5)
# yticks = np.linspace(0, ymax, 5)
# plt.figure(figsize=(4*cm, 4*cm))
# vec_logfc = df_protein['logFC'].values
# vec_log_pval = df_protein['LogPVal'].values
# mask = vec_log_pval > -np.log10(thresh_pval)
# mask_pos = mask * (vec_logfc > 0)
# mask_neg = mask * (vec_logfc < 0)
# gene_pos = df_protein['gene'].values[mask_pos]
# gene_neg = df_protein['gene'].values[mask_neg]
# plt.scatter(vec_logfc[~mask], vec_log_pval[~mask], s=0.5, c='gray', alpha=0.5, edgecolors='none', rasterized=True)
# plt.scatter(vec_logfc[mask_pos], vec_log_pval[mask_pos], s=0.5, c='red', edgecolors='none', rasterized=True)
# plt.scatter(vec_logfc[mask_neg], vec_log_pval[mask_neg], s=0.5, c='blue', edgecolors='none', rasterized=True)
# plt.text(xmax*0.9, ymax*0.9, '\n'.join(gene_pos[:10]), ha='right', va='top', c='red')
# plt.text(-xmax*0.9, ymax*0.9, '\n'.join(gene_neg[:10]), ha='left', va='top', c='blue')
# # plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]) + f'\n...\nTotal: {len(gene_pos)}', ha='right', va='top', c='red')
# # plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]) + f'\n...\nTotal: {len(gene_neg)}', ha='left', va='top', c='blue')
# # plt.title(f'{len(gene_pos)} up, {len(gene_neg)} down')
# plt.xlim([-xmax, xmax])
# plt.ylim([0, 1.05])
# plt.xticks(xticks)
# plt.yticks(yticks)
# # plt.xlabel('log2fc gene expression')
# # plt.ylabel('1 - FDR')
# plt.savefig(os.path.join(img_out, f'volcano_protein_{ds}_FDR{thresh_fdr}.{FMT}'), **fig_kwargs)

# df_protein = df_protein[df_protein['LogPVal'] > -np.log10(thresh_pval)]
# gene_protein_m6A_psi = []
# for this_gene, this_protein in df_protein[['gene', 'logFC']].values:
#     sub_df_logit = df_logit[df_logit['gene'] == this_gene]
#     if 'm6A' in sub_df_logit['mod'].values:
#         this_m6A = sub_df_logit[sub_df_logit['mod'] == 'm6A']['delta_logit'].values[0]
#     else:
#         this_m6A = np.nan
#     if 'psi' in sub_df_logit['mod'].values:
#         this_psi = sub_df_logit[sub_df_logit['mod'] == 'psi']['delta_logit'].values[0]
#     else:
#         this_psi = np.nan
#     gene_protein_m6A_psi.append(
#         (
#             this_gene,
#             this_protein,
#             this_m6A,
#             this_psi
#         )
#     )
# vec_gene, vec_protein, vec_m6A, vec_psi = np.vstack(gene_protein_m6A_psi).T
# vec_protein = np.float64(vec_protein)
# vec_m6A = np.float64(vec_m6A)
# vec_psi = np.float64(vec_psi)
#
# thresh_protein = 0
#
# plt.figure(figsize=(8*cm, 8*cm))
#
# plt.subplot(2, 2, 1)
# norm_pos_m6A, norm_neg_m6A, centers, width = plot_distribution(vec_m6A, vec_protein, thresh_protein,
#                                                                'm6A', bin_max=0.5, ymax=0.8)
# plt.subplot(2, 2, 2)
# plot_prob_diff(norm_pos_m6A, norm_neg_m6A, 'm6A', centers, width, ymax=0.4)
#
# plt.subplot(2, 2, 3)
# norm_pos_psi, norm_neg_psi, centers, width = plot_distribution(vec_psi, vec_protein, thresh_protein,
#                                                                'psi', bin_max=0.5, ymax=0.8)
# plt.subplot(2, 2, 4)
# plot_prob_diff(norm_pos_psi, norm_neg_psi, 'psi', centers, width, ymax=0.4)
#
# plt.savefig(os.path.join(img_out, f'hist_protein_vs_delta_logit_S_{ds}_pval{thresh_pval}.{FMT}'),
#             **fig_kwargs)