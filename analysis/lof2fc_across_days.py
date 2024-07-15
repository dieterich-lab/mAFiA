import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
mod_color = {
    'm6A': 'r',
    'psi': 'b'
}

res_dir = '/home/adrian/img_out/log2fc_per_gene'

days = ['day1', 'day7', 'day21', 'day56']
dfs = {}
for this_day in days:
     df_day = pd.read_csv(os.path.join(res_dir, f'gene_log2fc_{this_day}.tsv'), sep='\t')
     df_day.rename(columns={
         'num_mod_sites': f'num_mod_sites_{this_day}',
         'SHAM': f'SHAM_{this_day}',
         'TAC': f'TAC_{this_day}',
         'log2fc_mod': f'log2fc_mod_{this_day}',
         'log2fc_protein': f'log2fc_protein_{this_day}'
     }, inplace=True)
     dfs[this_day] = df_day

df_merged = pd.merge(dfs['day1'], dfs['day7'], on=['gene', 'mod'])
df_merged = pd.merge(df_merged, dfs['day21'], on=['gene', 'mod'])
df_merged = pd.merge(df_merged, dfs['day56'], on=['gene', 'mod'])

unique_genes = df_merged['gene'].unique()

ts_days = np.arange(1, 1+len(days))

gene_corr_mod_protein = {}
# this_gene = 'Fhl1'
for this_gene in unique_genes:
    df_gene = df_merged[df_merged['gene']==this_gene]
    if 'm6A' in df_gene['mod'].values:
        ts_log2fc_m6A = df_gene[df_gene['mod']=='m6A'][['log2fc_mod_day1', 'log2fc_mod_day7', 'log2fc_mod_day21', 'log2fc_mod_day56']].values[0]
    else:
        ts_log2fc_m6A = None
    if 'psi' in df_gene['mod'].values:
        ts_log2fc_psi = df_gene[df_gene['mod']=='psi'][['log2fc_mod_day1', 'log2fc_mod_day7', 'log2fc_mod_day21', 'log2fc_mod_day56']].values[0]
    else:
        ts_log2fc_psi = None
    ts_log2fc_protein = df_gene[['log2fc_protein_day1', 'log2fc_protein_day7', 'log2fc_protein_day21', 'log2fc_protein_day56']].values[0]

    if ts_log2fc_m6A is not None:
        corr_m6A_protein = np.corrcoef(ts_log2fc_m6A, ts_log2fc_protein)[0, 1]
    else:
        corr_m6A_protein = None
    if ts_log2fc_psi is not None:
        corr_psi_protein = np.corrcoef(ts_log2fc_psi, ts_log2fc_protein)[0, 1]
    else:
        corr_psi_protein = None
    gene_corr_mod_protein[this_gene] = (corr_m6A_protein, corr_psi_protein)

    if (ts_log2fc_m6A is not None) and (ts_log2fc_psi is not None):
        corr_m6A_psi = np.corrcoef(ts_log2fc_m6A, ts_log2fc_psi)[0, 1]

    if (np.nanmax((np.abs(corr_m6A_protein), np.abs(corr_psi_protein)))>=0.5) and (np.abs(corr_m6A_psi)>=0.5):
        plt.figure(figsize=(5, 4))
        plt.plot(ts_days, ts_log2fc_m6A, 'r-o', label=f"${{{dict_mod_display['m6A']}}}$, {corr_m6A_protein:.2f}")
        plt.plot(ts_days, ts_log2fc_psi, 'b-o', label=f"${{{dict_mod_display['psi']}}}$, {corr_psi_protein:.2f}")
        plt.plot(ts_days, ts_log2fc_protein, 'g-o', label='protein')
        plt.xticks(ts_days, days)
        plt.ylabel(f'$log_{2}FC$', fontsize=12)
        plt.legend(fontsize=10)
        plt.title(this_gene, fontsize=15)
