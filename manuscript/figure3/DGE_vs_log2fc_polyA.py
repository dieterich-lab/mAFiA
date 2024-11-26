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
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
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
gene_count_dir = os.path.join(res_dir, 'gene_counts')
polyA_dir = os.path.join(res_dir, 'polyA')

mods = ['m6A', 'psi']

thresh_pval = 0.05

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

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
df_abundance = df_abundance[df_abundance['PValue'] < thresh_pval]


### polyA ###
df_polyA = pd.read_csv(
    os.path.join(polyA_dir, f'gene_polyA_log2fc_pval_{ds}.tsv'),
    sep='\t'
)
df_polyA.rename(columns={'log2fc': 'log2fc_polyA'}, inplace=True)
df_polyA = df_polyA[df_polyA['pval'] < thresh_pval]

common_genes = set(df_abundance['gene']).intersection(set(df_polyA['gene']))
gene_dge_polyA = {}
for this_gene in common_genes:
    this_dge = df_abundance[df_abundance['gene'] == this_gene]['log2fc_gene'].values[0]
    this_polyA = df_polyA[df_polyA['gene'] == this_gene]['log2fc_polyA'].values[0]
    gene_dge_polyA[this_gene] = (this_dge, this_polyA)

vec_dge, vec_polyA = np.vstack(gene_dge_polyA.values()).T

# plt.figure()
# plt.scatter(vec_polyA, vec_dge, s=2)
# plt.xlim([-2, 2])
# plt.ylim([-2, 2])

ylim = [-0.5, 0.5]
yticks = np.round(np.linspace(*ylim, 3), 2)
plt.figure(figsize=(4*cm, 4*cm))
binned_y = [
    [x for x in vec_dge[vec_polyA < 0] if ~np.isnan(x) and ~np.isinf(x)],
    [x for x in vec_dge[vec_polyA > 0] if ~np.isnan(x) and ~np.isinf(x)],
]
num_genes = [len(this_bin) for this_bin in binned_y]
plt.boxplot(binned_y,
            positions=[-0.5, 0.5],
            widths=0.25,
            whis=0.5,
            showfliers=False)
plt.xticks([-0.5, 0.5],
           [rf'$\downarrow$ ({len(binned_y[0])})', rf'$\uparrow$ ({len(binned_y[1])})'])
# plt.yticks(yticks)
plt.xlabel('Delta polyA')
plt.ylabel('DGE')
plt.title(ds)
plt.savefig(os.path.join(img_out, f'boxplot_DGE_vs_log2fc_polyA_{ds}_pval{thresh_pval}.{FMT}'),
            **fig_kwargs)