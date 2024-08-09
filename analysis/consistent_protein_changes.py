import os
from functools import reduce
import pandas as pd
import re
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def load_df_gene():
    gene_bed_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'gene_id',
        'score',
        'strand',
        'source',
        'biotype',
        'other',
        'description'
    ]
    gene_bed_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
    df_gene = pd.read_csv(gene_bed_file, sep='\t', names=gene_bed_fields)
    df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
    df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]
    return df_gene

regulation = 'up'
days = ['day1', 'day7', 'day21', 'day56']

res_dir = '/home/adrian/img_out/protein_abundance'
ks_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/KS_test'
dfs_ks = {
    this_day: pd.read_csv(os.path.join(ks_dir, f'SHAM_TAC_{this_day}.bed'), sep='\t').drop(columns=['score'])
    for this_day in days
}

day_df = {
    this_day: pd.read_csv(os.path.join(res_dir, f'protein_{regulation}_{this_day}.tsv'), sep='\t')
    for this_day in days
}

df_merged = reduce(lambda left, right: pd.merge(left, right, on=['gene'], how='outer'), list(day_df.values()))
df_summary = df_merged[['gene', 'percentage_change_day1', 'percentage_change_day7', 'percentage_change_day21', 'percentage_change_day56']]

df_gene_region = load_df_gene()

this_gene = 'Fhl1'
thresh_pval = 1.0

this_gene_df = df_gene_region[df_gene_region['gene'] == this_gene]
this_chrom = this_gene_df['chrom'].unique()[0]
this_chromStart = this_gene_df['chromStart'].min()
this_chromEnd = this_gene_df['chromEnd'].max()
this_strand = this_gene_df['strand'].unique()[0]
this_gene_dfs_ks = {}
for this_day in days:
    this_day_df_ks = dfs_ks[this_day]
    this_gene_dfs_ks[this_day] = this_day_df_ks[
        (this_day_df_ks['chrom'] == this_chrom)
        * (this_day_df_ks['chromStart'] >= this_chromStart)
        * (this_day_df_ks['chromEnd'] <= this_chromEnd)
        * (this_day_df_ks['strand'] == this_strand)
        # * (this_day_df_ks['pval'] < thresh_pval)
    ]

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

ylim = [-12, 12]
plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    day_deltas = []
    for this_day in days:
        this_day_df = this_gene_dfs_ks[this_day]
        if len(this_day_df):
            day_deltas.append(this_day_df[this_day_df['name'] == this_mod]['delta'].values)
        else:
            day_deltas.append([0])
    plt.violinplot(day_deltas, range(len(day_deltas)), showmedians=True)
    plt.xticks(range(len(day_deltas)), days)
    plt.ylabel(rf'$\Delta S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
    plt.ylim(ylim)
plt.suptitle(this_gene, fontsize=15)
