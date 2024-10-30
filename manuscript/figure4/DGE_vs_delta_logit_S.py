import os
from glob import glob
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
abundance_dir = os.path.join(res_dir, 'DTU')
logit_dir = os.path.join(res_dir, 'transcript_logit')

mods = ['m6A', 'psi']

thresh_fdr = 0.1

# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']
# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']
ds = 'Diet'
conditions = ['WT_CD_merged', 'WT_WD_merged']

df_abundance = pd.read_excel(
    glob(os.path.join(abundance_dir, f'DTU_{ds}_*.xlsx'))[0],
    'DGE'
)
df_abundance = df_abundance[df_abundance['FDR'] < thresh_fdr]
df_abundance.rename(columns={'mgi_symbol': 'gene'}, inplace=True)

df_logit = pd.read_csv(
    os.path.join(logit_dir, f'delta_logitS_{ds}.tsv'),
    sep='\t'
)
# df_logit = df_logit[df_logit['pval'] < 0.05]

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
bin_max = 1.5
num_bins = 10
bin_edges = np.round(np.linspace(-bin_max, bin_max, num_bins+1), 2)
bin_width = np.round(bin_edges[1] - bin_edges[0], 2)
bin_centers = np.round(0.5 * (bin_edges[1:] + bin_edges[:-1]), 2)

plt.figure(figsize=(10, 10))

plt.subplot(2, 2, 1)
hist_pos, _ = np.histogram(vec_m6A[vec_dge > thresh_dge], range=[-bin_max, bin_max], bins=num_bins)
num_pos = np.sum(hist_pos)
norm_pos = hist_pos / num_pos
hist_neg, _ = np.histogram(vec_m6A[vec_dge < -thresh_dge], range=[-bin_max, bin_max], bins=num_bins)
num_neg = np.sum(hist_neg)
norm_neg = hist_neg / num_neg
plt.bar(bin_centers, norm_pos, label=f'DGE positive ({num_pos})', width=bin_width*0.8, alpha=0.5)
plt.bar(bin_centers, norm_neg, label=f'DGE negative ({num_neg})', width=bin_width*0.8, alpha=0.5)
plt.legend(loc='upper left')
plt.xlabel('delta logit S')
plt.ylabel('density')
plt.title('m6A')

plt.subplot(2, 2, 2)
plt.bar(bin_centers, norm_pos-norm_neg, width=bin_width*0.8)
plt.axhline(y=0, c='gray', ls='--')
plt.xlabel('delta logit S')
plt.ylabel('pr(DGE pos) - pr(DGE neg)')
plt.title('m6A')

plt.subplot(2, 2, 3)
hist_pos, _ = np.histogram(vec_psi[vec_dge > thresh_dge], range=[-bin_max, bin_max], bins=num_bins)
num_pos = np.sum(hist_pos)
norm_pos = hist_pos / num_pos
hist_neg, _ = np.histogram(vec_psi[vec_dge < -thresh_dge], range=[-bin_max, bin_max], bins=num_bins)
num_neg = np.sum(hist_neg)
norm_neg = hist_neg / num_neg
plt.bar(bin_centers, norm_pos, label=f'DGE positive ({num_pos})', width=bin_width*0.8, alpha=0.5)
plt.bar(bin_centers, norm_neg, label=f'DGE negative ({num_neg})', width=bin_width*0.8, alpha=0.5)
plt.legend(loc='upper left')
plt.xlabel('delta logit S')
plt.ylabel('density')
plt.title('psi')

plt.subplot(2, 2, 4)
plt.bar(bin_centers, norm_pos-norm_neg, width=bin_width*0.8)
plt.axhline(y=0, c='gray', ls='--')
plt.xlabel('delta logit S')
plt.ylabel('pr(DGE pos) - pr(DGE neg)')
plt.title('psi')

plt.suptitle(ds)