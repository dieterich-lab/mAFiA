import os
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import pysam
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import re

dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

day = 'day21'
thresh_pval = 1.0
num_bases = 5

lower_limit = -0.5
upper_limit = 0.5


proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'
df_proteome = pd.read_csv(proteome_file, sep='\t')
df_proteome_sel = df_proteome[
    (df_proteome['tp']==f"Day {day.lstrip('day')}")
    * (df_proteome['p.adj']<thresh_pval)
    * (df_proteome['main.comparison']=='TAC vs Sham')
    ]

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
df_gene['gene_name'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

proteome_conditions = [
    f'logFC $< {str(lower_limit)}$',
    # 'logFC $\smallin [-0.5, 0.5)$',
    f'logFC $\geq {str(upper_limit)}$'
]
proteome_colors = [
    'skyblue',
    # 'royalblue',
    'midnightblue'
]
bam_files = {
    'TAC': f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_{day}/chrALL.mAFiA.reads.bam',
    'SHAM': f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_{day}/chrALL.mAFiA.reads.bam',
}

num_bins = 20
max_bin = 1.0

img_out = '/home/adrian/img_out/proteome_avg_mod_per_read'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def calc_mod_level_per_read_aligned_to_region(in_bam_file, in_chrom, in_chromStart, in_chromEnd, in_strand, num_bases):
    read_avg_mod_level = {}
    if in_strand=='+':
        flag_require = 0
    else:
        flag_require = 16
    with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
        for this_read in in_bam.fetch(in_chrom, in_chromStart, in_chromEnd):
            if this_read.flag == flag_require:
                read_avg_mod_level[this_read.query_name] = {}
                for this_mod in dict_mod_code.keys():
                    mod_bases = this_read.modified_bases.get(('N', 0, dict_mod_code[this_mod]), [])
                    if len(mod_bases)>=num_bases:
                        read_avg_mod_level[this_read.query_name][this_mod] = np.mean(np.partition([tup[1] for tup in mod_bases], -num_bases)[-num_bases:]) / 255.0
                    else:
                        read_avg_mod_level[this_read.query_name][this_mod] = np.mean([tup[1] for tup in mod_bases]) / 255.0
    return read_avg_mod_level


def collect_gene_avg_mod_level(in_df_proteome, in_bam_file):
    gene_avg_mod_level = {}
    for this_gene in tqdm(in_df_proteome['SYMBOL']):
        sub_df_gene = df_gene[df_gene['gene_name'] == this_gene]
        if len(sub_df_gene) == 0:
            continue
        this_chrom = sub_df_gene['chrom'].iloc[0]
        this_chromStart = sub_df_gene['chromStart'].min()
        this_chromEnd = sub_df_gene['chromEnd'].max()
        this_strand = sub_df_gene['strand'].iloc[0]
        gene_avg_mod_level[this_gene] = calc_mod_level_per_read_aligned_to_region(
            in_bam_file, this_chrom, this_chromStart, this_chromEnd, this_strand, num_bases=num_bases)
    return gene_avg_mod_level


df_proteome_negative = df_proteome_sel[df_proteome_sel['logFC'] < lower_limit]
# df_proteome_neutral = df_proteome_sel[(df_proteome_sel['logFC'] >= -0.5) * (df_proteome_sel['logFC'] < 0.5)]
df_proteome_positive = df_proteome_sel[df_proteome_sel['logFC'] >= upper_limit]
cond_df_proteome = {
    this_cond: this_df for this_cond, this_df in
    zip(proteome_conditions,
        [df_proteome_negative,
         # df_proteome_neutral,
         df_proteome_positive])
}

experiments = ['SHAM', 'TAC']
exp_ls = {
    'SHAM': '--',
    'TAC': '-'
}

plt.figure(figsize=(10, 10))
for row_ind, this_exp in enumerate(experiments):
    for col_ind, this_mod in enumerate(mods):
        for this_proteome_cond, this_color in zip(proteome_conditions, proteome_colors):
            plt.subplot(2, 2, row_ind*2+col_ind+1)
            all_gene_avg_mod_level = collect_gene_avg_mod_level(cond_df_proteome[this_proteome_cond], bam_files[this_exp])
            this_mod_vals = [this_read_mod_level.get(this_mod, np.nan)
                             for this_gene_read_mod_level in all_gene_avg_mod_level.values()
                             for this_read_mod_level in this_gene_read_mod_level.values()]
            this_hist, bin_edges = np.histogram(this_mod_vals, range=[0, max_bin], bins=num_bins)
            this_hist = this_hist / np.sum(this_hist)
            bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
            plt.plot(bin_centers, this_hist, c=this_color, label=this_proteome_cond)
        plt.xlim([-0.01, 1.01])
        # plt.ylim([0, 0.1])
        plt.legend(fontsize=8)
        plt.title(this_exp, fontsize=12)
        plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $ per read', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
plt.suptitle(day, fontsize=15)
plt.savefig(os.path.join(img_out, f'hist_avg_mods_per_read_grouped_by_protein_expression_{day}.png'),
            bbox_inches='tight')
