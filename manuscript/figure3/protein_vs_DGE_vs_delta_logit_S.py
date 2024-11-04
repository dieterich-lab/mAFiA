import os
from glob import glob
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import numpy as np
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


img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
abundance_dir = os.path.join(res_dir, 'DTU')
logit_dir = os.path.join(res_dir, 'transcript_logit')
gene_count_dir = os.path.join(res_dir, 'gene_counts')
protein_dir = os.path.join(res_dir, 'protein_abundance')

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

### volcano plot ###
xmax = 4
xticks = np.linspace(-xmax, xmax, 5)
yticks = np.linspace(0, 1, 5)
plt.figure(figsize=(4*cm, 4*cm))
vec_logfc = df_abundance['logFC'].values
vec_inv_fdr = 1.0 - df_abundance['FDR'].values
mask = vec_inv_fdr > (1.0 - thresh_fdr)
mask_pos = mask * (vec_logfc > 0)
mask_neg = mask * (vec_logfc < 0)
gene_pos = df_abundance['gene'].values[mask_pos]
gene_neg = df_abundance['gene'].values[mask_neg]
plt.scatter(vec_logfc[~mask], vec_inv_fdr[~mask], s=0.5, c='gray', alpha=0.5, edgecolors='none', rasterized=True)
plt.scatter(vec_logfc[mask_pos], vec_inv_fdr[mask_pos], s=0.5, c='red', edgecolors='none', rasterized=True)
plt.scatter(vec_logfc[mask_neg], vec_inv_fdr[mask_neg], s=0.5, c='blue', edgecolors='none', rasterized=True)
plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]), ha='right', va='top', c='red')
plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]), ha='left', va='top', c='blue')
# plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]) + f'\n...\nTotal: {len(gene_pos)}', ha='right', va='top', c='red')
# plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]) + f'\n...\nTotal: {len(gene_neg)}', ha='left', va='top', c='blue')
# plt.title(f'{len(gene_pos)} up, {len(gene_neg)} down')
plt.xlim([-xmax, xmax])
plt.ylim([0, 1.05])
plt.xticks(xticks)
plt.yticks(yticks)
# plt.xlabel('log2fc gene expression')
# plt.ylabel('1 - FDR')
plt.savefig(os.path.join(img_out, f'volcano_DGE_{ds}_FDR{thresh_fdr}.{FMT}'), **fig_kwargs)

df_abundance = df_abundance[df_abundance['FDR'] < thresh_fdr]
df_abundance = df_abundance[['BambuGene' not in this_id for this_id in df_abundance['ID']]]

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
out_abundance_filename = os.path.join(abundance_dir, f'DGE_{ds}_with_gene_count.tsv')
df_abundance.to_csv(out_abundance_filename, sep='\t', index=False)

thresh_coverage = 20
num_sufficient_coverage = (df_abundance[
                               [f'gene_count_{this_cond}' for this_cond in conditions]
                           ].values >= thresh_coverage).all(axis=1).sum()
with open(out_abundance_filename, 'a') as fout:
    fout.write(f'\n\n#{ds}: {num_sufficient_coverage} / {len(df_abundance)} with coverage >= {thresh_coverage}')

df_logit = pd.read_csv(
    os.path.join(logit_dir, f'delta_logitS_{ds}.tsv'),
    sep='\t'
)

### volcano plot ###
xmax = 4
ymax = 20
xticks = np.linspace(-xmax, xmax, 5)
yticks = np.linspace(0, ymax, 5)

plt.figure(figsize=(8*cm, 4*cm))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    df_mod = df_logit[df_logit['mod'] == this_mod]
    df_mod.sort_values('pval', inplace=True)
    vec_delta_logit = df_mod['delta_logit'].values
    vec_neg_log_pval = -np.log10(df_mod['pval'].values)
    mask = vec_neg_log_pval > -np.log10(thresh_pval)
    mask_pos = mask * (vec_delta_logit > 0)
    mask_neg = mask * (vec_delta_logit < 0)
    gene_pos = df_mod['gene'].values[mask_pos]
    gene_neg = df_mod['gene'].values[mask_neg]
    plt.scatter(vec_delta_logit[~mask], vec_neg_log_pval[~mask], s=0.5, c='gray', alpha=0.5, edgecolors='none', rasterized=True)
    plt.scatter(vec_delta_logit[mask_pos], vec_neg_log_pval[mask_pos], s=0.5, c='red', edgecolors='none', rasterized=True)
    plt.scatter(vec_delta_logit[mask_neg], vec_neg_log_pval[mask_neg], s=0.5, c='blue', edgecolors='none', rasterized=True)
    plt.text(xmax*0.9, ymax*0.9, '\n'.join(gene_pos[:10]), ha='right', va='top', c='red')
    plt.text(-xmax*0.9, ymax*0.9, '\n'.join(gene_neg[:10]), ha='left', va='top', c='blue')
    # plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]) + f'\n...\nTotal: {len(gene_pos)}', ha='right', va='top', c='red')
    # plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]) + f'\n...\nTotal: {len(gene_neg)}', ha='left', va='top', c='blue')
    # plt.title(f'{len(gene_pos)} up, {len(gene_neg)} down')
    plt.xlim([-xmax, xmax])
    plt.ylim([0, ymax])
    plt.xticks(xticks)
    plt.yticks(yticks)
# plt.xlabel('log2fc gene expression')
# plt.ylabel('1 - FDR')
plt.savefig(os.path.join(img_out, f'volcano_delta_logit_{ds}_pval{thresh_pval}.{FMT}'), **fig_kwargs)

df_logit = df_logit[df_logit['pval'] < thresh_pval]

gene_dge_m6A_psi = []
for this_gene, this_dge in df_abundance[['gene', 'logFC']].values:
    sub_df_logit = df_logit[df_logit['gene'] == this_gene]
    if 'm6A' in sub_df_logit['mod'].values:
        this_m6A = sub_df_logit[sub_df_logit['mod'] == 'm6A']['delta_logit'].values[0]
    else:
        this_m6A = np.nan
    if 'psi' in sub_df_logit['mod'].values:
        this_psi = sub_df_logit[sub_df_logit['mod'] == 'psi']['delta_logit'].values[0]
    else:
        this_psi = np.nan
    gene_dge_m6A_psi.append(
        (
            this_gene,
            this_dge,
            this_m6A,
            this_psi
        )
    )

vec_gene, vec_dge, vec_m6A, vec_psi = np.vstack(gene_dge_m6A_psi).T
vec_dge = np.float64(vec_dge)
vec_m6A = np.float64(vec_m6A)
vec_psi = np.float64(vec_psi)

thresh_dge = 0

plt.figure(figsize=(8*cm, 8*cm))
plt.subplot(2, 2, 1)
norm_pos_m6A, norm_neg_m6A, centers, width = plot_distribution(vec_m6A, vec_dge, thresh_dge, 'm6A',
                                                               bin_max=0.5, ymax=0.7)
plt.subplot(2, 2, 2)
plot_prob_diff(norm_pos_m6A, norm_neg_m6A, 'm6A', centers, width, ymax=0.6)

plt.subplot(2, 2, 3)
norm_pos_psi, norm_neg_psi, centers, width = plot_distribution(vec_psi, vec_dge, thresh_dge, 'psi',
                                                               bin_max=0.5, ymax=0.7)
plt.subplot(2, 2, 4)
plot_prob_diff(norm_pos_psi, norm_neg_psi, 'psi', centers, width, ymax=0.6)

# plt.suptitle(ds)
plt.savefig(os.path.join(img_out, f'hist_DGE_vs_delta_logit_S_{ds}_FDR{thresh_fdr}_pval{thresh_pval}.{FMT}'),
            **fig_kwargs)

### protein abundance ###
df_protein = pd.read_csv(
    os.path.join(protein_dir, 'TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'),
    sep='\t'
)
df_protein = df_protein[df_protein['tp'] == 'Combined']
df_protein.sort_values('p.adj', inplace=True)
df_protein.rename(columns={'SYMBOL': 'gene'}, inplace=True)

### volcano plot ###
xmax = 2
ymax = 8
xticks = np.linspace(-xmax, xmax, 5)
yticks = np.linspace(0, ymax, 5)
plt.figure(figsize=(4*cm, 4*cm))
vec_logfc = df_protein['logFC'].values
vec_log_pval = df_protein['LogPVal'].values
mask = vec_log_pval > -np.log10(thresh_pval)
mask_pos = mask * (vec_logfc > 0)
mask_neg = mask * (vec_logfc < 0)
gene_pos = df_protein['gene'].values[mask_pos]
gene_neg = df_protein['gene'].values[mask_neg]
plt.scatter(vec_logfc[~mask], vec_log_pval[~mask], s=0.5, c='gray', alpha=0.5, edgecolors='none', rasterized=True)
plt.scatter(vec_logfc[mask_pos], vec_log_pval[mask_pos], s=0.5, c='red', edgecolors='none', rasterized=True)
plt.scatter(vec_logfc[mask_neg], vec_log_pval[mask_neg], s=0.5, c='blue', edgecolors='none', rasterized=True)
plt.text(xmax*0.9, ymax*0.9, '\n'.join(gene_pos[:10]), ha='right', va='top', c='red')
plt.text(-xmax*0.9, ymax*0.9, '\n'.join(gene_neg[:10]), ha='left', va='top', c='blue')
# plt.text(xmax * 0.9, 0.9, '\n'.join(gene_pos[:10]) + f'\n...\nTotal: {len(gene_pos)}', ha='right', va='top', c='red')
# plt.text(-xmax * 0.9, 0.9, '\n'.join(gene_neg[:10]) + f'\n...\nTotal: {len(gene_neg)}', ha='left', va='top', c='blue')
# plt.title(f'{len(gene_pos)} up, {len(gene_neg)} down')
plt.xlim([-xmax, xmax])
plt.ylim([0, 1.05])
plt.xticks(xticks)
plt.yticks(yticks)
# plt.xlabel('log2fc gene expression')
# plt.ylabel('1 - FDR')
plt.savefig(os.path.join(img_out, f'volcano_protein_{ds}_FDR{thresh_fdr}.{FMT}'), **fig_kwargs)

df_protein = df_protein[df_protein['LogPVal'] > -np.log10(thresh_pval)]
gene_protein_m6A_psi = []
for this_gene, this_protein in df_protein[['gene', 'logFC']].values:
    sub_df_logit = df_logit[df_logit['gene'] == this_gene]
    if 'm6A' in sub_df_logit['mod'].values:
        this_m6A = sub_df_logit[sub_df_logit['mod'] == 'm6A']['delta_logit'].values[0]
    else:
        this_m6A = np.nan
    if 'psi' in sub_df_logit['mod'].values:
        this_psi = sub_df_logit[sub_df_logit['mod'] == 'psi']['delta_logit'].values[0]
    else:
        this_psi = np.nan
    gene_protein_m6A_psi.append(
        (
            this_gene,
            this_protein,
            this_m6A,
            this_psi
        )
    )
vec_gene, vec_protein, vec_m6A, vec_psi = np.vstack(gene_protein_m6A_psi).T
vec_protein = np.float64(vec_protein)
vec_m6A = np.float64(vec_m6A)
vec_psi = np.float64(vec_psi)

thresh_protein = 0

plt.figure(figsize=(8*cm, 8*cm))

plt.subplot(2, 2, 1)
norm_pos_m6A, norm_neg_m6A, centers, width = plot_distribution(vec_m6A, vec_protein, thresh_protein,
                                                               'm6A', bin_max=0.5, ymax=0.8)
plt.subplot(2, 2, 2)
plot_prob_diff(norm_pos_m6A, norm_neg_m6A, 'm6A', centers, width, ymax=0.4)

plt.subplot(2, 2, 3)
norm_pos_psi, norm_neg_psi, centers, width = plot_distribution(vec_psi, vec_protein, thresh_protein,
                                                               'psi', bin_max=0.5, ymax=0.8)
plt.subplot(2, 2, 4)
plot_prob_diff(norm_pos_psi, norm_neg_psi, 'psi', centers, width, ymax=0.4)

plt.savefig(os.path.join(img_out, f'hist_protein_vs_delta_logit_S_{ds}_pval{thresh_pval}.{FMT}'),
            **fig_kwargs)