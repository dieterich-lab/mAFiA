import os
import ast
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import re
import numpy as np
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

ds = 'HFpEF'
conditions = ['ctrl', 'HFpEF']
cond_colors = {
    'ctrl': 'b',
    'HFpEF': 'r'
}

gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'frame',
    'description'
]
df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

cond_mod_profile = {this_cond: {} for this_cond in conditions}
# this_cond = 'ctrl'

num_bins = 100
bin_edges = np.linspace(0, 1, num_bins+1)
bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

thresh_conf = 0.0

for this_cond in conditions:
    df_sites = pd.read_csv(
        os.path.join(res_dir, ds, f'{this_cond}_merged', 'chrALL.mAFiA.sites.bed'),
        sep='\t', dtype={'chrom': str}
    )
    df_sites_thresh = df_sites[df_sites['confidence'] >= thresh_conf]

    mod_gene_normPos_modRatio = {this_mod: {} for this_mod in mods}
    for this_gene in tqdm(df_gene['gene']):
        sub_df_gene = df_gene[df_gene['gene'] == this_gene]
        chrom = sub_df_gene['chrom'].unique()[0]
        geneStart, geneEnd = sub_df_gene[['chromStart', 'chromEnd']].values[0]
        gene_len = geneEnd - geneStart
        strand = sub_df_gene['strand'].unique()[0]

        sub_df_sites = df_sites_thresh[
            (df_sites_thresh['chrom']==chrom)
            * (df_sites_thresh['chromStart'] >= geneStart)
            * (df_sites_thresh['chromEnd'] <= geneEnd)
            * (df_sites_thresh['strand'] == strand)
        ]

        if len(sub_df_sites):
            for this_mod in mods:
                mod_sub_df_sites = sub_df_sites[sub_df_sites['name'] == this_mod]
                if len(mod_sub_df_sites):
                    vec_refPos = mod_sub_df_sites['chromStart'].values
                    vec_modRatio = mod_sub_df_sites['modRatio'].values
                    if strand == '+':
                        vec_normPos = (vec_refPos - geneStart) / gene_len
                    else:
                        vec_normPos = (geneEnd - vec_refPos - 1) / gene_len

                    mod_gene_normPos_modRatio[this_mod][this_gene] = list(zip(vec_normPos, vec_modRatio))

#     for this_mod in mods:
#         this_mod_hist, _ = np.histogram(thresh_mod_pos[this_mod], range=[0, 1], bins=num_bins)
#         cond_mod_profile[this_cond][this_mod] = this_mod_hist
#
# plt.figure(figsize=(10, 4))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     for this_cond in conditions:
#         this_hist = cond_mod_profile[this_cond][this_mod]
#         this_hist = this_hist / np.sum(this_hist)
#         plt.plot(bin_centers, this_hist, c=cond_colors[this_cond], label=this_cond)
#     plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
