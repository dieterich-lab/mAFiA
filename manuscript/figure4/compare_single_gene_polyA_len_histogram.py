import os
from glob import glob
import pandas as pd
from collections import Counter
import numpy as np
from scipy.stats import kstest
from tqdm import tqdm
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


def get_df_log2fc(in_cond_df):
    common_genes = list(set(in_cond_df[conditions[0]]['gene'].unique()).intersection(set(in_cond_df[conditions[1]]['gene'].unique())))
    out_gene_polyA_log2fc_pval = {}
    for thisGene in tqdm(common_genes):
        polyA_len_0 = in_cond_df[conditions[0]][in_cond_df[conditions[0]]['gene'] == thisGene]['polyA_length']
        polyA_len_1 = in_cond_df[conditions[1]][in_cond_df[conditions[1]]['gene'] == thisGene]['polyA_length']
        mean_len_0 = np.mean(polyA_len_0)
        mean_len_1 = np.mean(polyA_len_1)
        log2fc = np.log2(mean_len_1 / mean_len_0)
        ks_stat, pval = kstest(polyA_len_0, polyA_len_1)
        out_gene_polyA_log2fc_pval[thisGene] = (len(polyA_len_0), len(polyA_len_1), mean_len_0, mean_len_1, log2fc, pval)

    out_df_gene_polyA_log2fc_pval = pd.DataFrame([[k] + list(v) for k, v in out_gene_polyA_log2fc_pval.items()],
                                             columns=['gene',
                                                      f'num_reads_{conditions[0]}',
                                                      f'num_reads_{conditions[1]}',
                                                      f'mean_polyA_len_{conditions[0]}',
                                                      f'mean_polyA_len_{conditions[1]}',
                                                      'log2fc',
                                                      'pval'
                                                      ])
    out_df_gene_polyA_log2fc_pval.sort_values('log2fc', inplace=True)
    return out_df_gene_polyA_log2fc_pval


polyA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'

ds = 'TAC'
conditions = ['SHAM', 'TAC']
# ds = 'HFpEF'
# conditions = ['ctrl', 'HFpEF']
# ds = 'Diet'
# conditions = ['WT_CD', 'WT_WD']

cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['b', 'r'])}

log2fc_filename = os.path.join(polyA_dir, f'gene_polyA_log2fc_pval_{ds}.tsv')

cond_df = {}
for this_cond in conditions:
    all_filepaths = glob(os.path.join(polyA_dir, f'polyA_reads_annotated_{ds}_{this_cond}_day*.tsv'))
    # all_filepaths = glob(os.path.join(polyA_dir, f'polyA_reads_annotated_{ds}_{this_cond}.tsv'))
    all_dfs = []
    for this_filepath in all_filepaths:
        all_dfs.append(pd.read_csv(this_filepath, sep='\t'))
    cond_df[this_cond] = pd.concat(all_dfs).sort_values('gene')

if os.path.exists(log2fc_filename):
    df_gene_polyA_log2fc_pval = pd.read_csv(log2fc_filename, sep='\t')
else:
    df_gene_polyA_log2fc_pval = get_df_log2fc(cond_df)
    df_gene_polyA_log2fc_pval.to_csv(log2fc_filename, sep='\t', index=False, float_format='%.6f')


### volcano plot ###
df_gene_polyA_log2fc_pval['logit'] = -np.log10(df_gene_polyA_log2fc_pval['pval'])

thresh_logit = 3
thresh_log2fc = 0.0

positive_mask = (df_gene_polyA_log2fc_pval['log2fc'].values >= thresh_log2fc) \
                * (df_gene_polyA_log2fc_pval['logit'].values >= thresh_logit)
negative_mask = (df_gene_polyA_log2fc_pval['log2fc'].values < -thresh_log2fc) \
                * (df_gene_polyA_log2fc_pval['logit'].values >= thresh_logit)
positive_genes = df_gene_polyA_log2fc_pval['gene'].values[positive_mask]
negative_genes = df_gene_polyA_log2fc_pval['gene'].values[negative_mask]

xmax = 2
ymax = 5.5
xticks = np.linspace(-xmax, xmax, 5)
yticks = np.arange(0, ymax+0.01, 1)

plt.figure(figsize=(4*cm, 4*cm))
plt.scatter(df_gene_polyA_log2fc_pval['log2fc'], df_gene_polyA_log2fc_pval['logit'], c='gray', s=1, alpha=0.5)
plt.scatter(df_gene_polyA_log2fc_pval['log2fc'].values[positive_mask],
            df_gene_polyA_log2fc_pval['logit'].values[positive_mask], c='red', s=1)
plt.scatter(df_gene_polyA_log2fc_pval['log2fc'].values[negative_mask],
            df_gene_polyA_log2fc_pval['logit'].values[negative_mask], c='blue', s=1)
plt.xlim([-2, 2])
plt.ylim([0, 5.5])
plt.xticks(xticks)
plt.yticks(yticks)
plt.text(-xmax+0.05, ymax-0.05, '\n'.join(negative_genes), c='b', ha='left', va='top')
plt.text(xmax-0.05, ymax-0.05, '\n'.join(positive_genes), c='r', ha='right', va='top')

plt.savefig(os.path.join(img_out, f'volcano_plot_polyA_len_log2fc_pval_{ds}.{FMT}'), **fig_kwargs)


### single transcript distribution ###
this_gene = 'Nppb'
gene_cond_distribution = {}
for this_cond in conditions:
    this_cond_df = cond_df[this_cond]
    gene_cond_distribution[this_cond] = this_cond_df[this_cond_df['gene'] == this_gene]['polyA_length']

xmax = 300
ymax = 0.4
xticks = np.linspace(0, xmax, 4)
yticks = np.arange(0, ymax+0.01, 0.1)
bin_width = 20
num_bins = np.int64(xmax / bin_width)
bin_edges = np.linspace(0, xmax, num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

plt.figure(figsize=(4*cm, 4*cm))
for this_cond in conditions:
    this_hist, _ = np.histogram(gene_cond_distribution[this_cond], bins=bin_edges)
    norm_hist = this_hist / np.sum(this_hist)
    plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=this_cond)
plt.xticks(xticks)
plt.yticks(yticks)
plt.xlim([0, xmax])
plt.ylim([0, ymax])
plt.savefig(os.path.join(img_out, f'hist_polyA_len_{ds}_{this_gene}.{FMT}'), **fig_kwargs)
