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

condition = 'TAC'
day = 'day7'
num_bases = 5

proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/protein_expression_data_mapped.tsv'
df_proteome = pd.read_csv(proteome_file, sep='\t')
sel_cols = list(df_proteome.columns[df_proteome.columns.str.contains(f'{condition}_{day}')])
df_proteome_sel = df_proteome[['gene'] + sel_cols]
df_proteome_sel[f'{condition}_{day}_avg'] = df_proteome_sel[sel_cols].mean(axis=1)

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

# proteome_conditions = [
#     f'$< {str(lower_limit)}$',
#     f'$\smallin [{lower_limit}, {upper_limit})$',
#     f'$\geq {str(upper_limit)}$'
# ]
proteome_conditions = [
    # 'q[0,50)',
    'q[0,25)',
    # 'q[25,50)',
    # 'q[50,75)',
    'q[25,75)',
    'q[75,100)'
    # 'q[50,100)'
]
proteome_colors = [
    'skyblue',
    'royalblue',
    'midnightblue',
    # 'darkblue'
]
bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/{condition}_{day}/chrALL.mAFiA.reads.bam'
polyA_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA/{condition}_{day}_read2polyA_length.txt'

num_bins = 20
max_bin = 1.0

img_out = '/home/adrian/img_out/proteome_avg_mod_per_read'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}


def collect_reads_aligned_to_region(in_bam_file, in_chrom, in_chromStart, in_chromEnd, in_strand):
    out_reads = []
    if in_strand=='+':
        flag_require = 0
    else:
        flag_require = 16
    with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
        for this_read in in_bam.fetch(in_chrom, in_chromStart, in_chromEnd):
            if this_read.flag == flag_require:
                out_reads.append(this_read)
    return out_reads


def collect_proteome_reads(in_df_proteome, in_bam_file):
    gene_reads = {}
    for this_gene in tqdm(in_df_proteome['gene']):
        sub_df_gene = df_gene[df_gene['gene_name'] == this_gene]
        if len(sub_df_gene) == 0:
            continue
        this_chrom = sub_df_gene['chrom'].iloc[0]
        this_chromStart = sub_df_gene['chromStart'].min()
        this_chromEnd = sub_df_gene['chromEnd'].max()
        this_strand = sub_df_gene['strand'].iloc[0]
        gene_reads[this_gene] = collect_reads_aligned_to_region(
            in_bam_file, this_chrom, this_chromStart, this_chromEnd, this_strand)
    return gene_reads


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


def calc_mod_level_per_read_in_all_genes(in_gene_reads, in_num_bases):
    gene_avg_mod_level = {}
    for this_gene, this_gene_reads in tqdm(in_gene_reads.items()):
        read_avg_mod_level = {}
        for this_read in this_gene_reads:
            read_avg_mod_level[this_read.query_name] = {}
            for mod in dict_mod_code.keys():
                mod_bases = this_read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])
                if len(mod_bases) >= in_num_bases:
                    read_avg_mod_level[this_read.query_name][mod] = np.mean(np.partition([tup[1] for tup in mod_bases], -in_num_bases)[-in_num_bases:]) / 255.0
                else:
                    read_avg_mod_level[this_read.query_name][this_mod] = np.mean([tup[1] for tup in mod_bases]) / 255.0
        if len(read_avg_mod_level):
            gene_avg_mod_level[this_gene] = read_avg_mod_level
    return gene_avg_mod_level


def collect_gene_avg_mod_level(in_df_proteome, in_bam_file):
    gene_avg_mod_level = {}
    for this_gene in tqdm(in_df_proteome['gene']):
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


def collect_polyA_lengths_from_read_ids(in_df_polyA, in_read_ids):
    return in_df_polyA[in_df_polyA['read_id'].isin(in_read_ids)]['polyA_len'].values

# df_proteome_negative = df_proteome_sel[df_proteome_sel[f'{condition}_{day}_avg'] < lower_limit]
# df_proteome_neutral = df_proteome_sel[(df_proteome_sel[f'{condition}_{day}_avg'] >= lower_limit)
#                                       * (df_proteome_sel[f'{condition}_{day}_avg'] < upper_limit)]
# df_proteome_positive = df_proteome_sel[df_proteome_sel[f'{condition}_{day}_avg'] >= upper_limit]
# cond_df_proteome = {
#     this_cond: this_df for this_cond, this_df in
#     zip(proteome_conditions,
#         [df_proteome_negative,
#          df_proteome_neutral,
#          df_proteome_positive])
# }

df_polyA = pd.read_csv(polyA_file, sep='\t', names=['read_id', 'polyA_len'])

cond_df_proteome = {}
for this_cond in proteome_conditions:
    q_low, q_high = this_cond.lstrip('q[').rstrip(')').split(',')
    p_low = df_proteome_sel[f'{condition}_{day}_avg'].quantile(int(q_low)/100.0)
    p_high = df_proteome_sel[f'{condition}_{day}_avg'].quantile(int(q_high)/100.0)
    if q_high=='100':
        p_high += 0.000001
    cond_df_proteome[this_cond] = df_proteome_sel[(df_proteome_sel[f'{condition}_{day}_avg'] >= p_low)
                                                  * (df_proteome_sel[f'{condition}_{day}_avg'] < p_high)]

cond_reads = {this_cond: collect_proteome_reads(this_df, bam_file) for this_cond, this_df in cond_df_proteome.items()}

plt.figure(figsize=(9, 16))

num_rows = 4
num_cols = 2

plt.subplot(num_rows, num_cols, 1)
for this_proteome_cond, this_color in zip(proteome_conditions, proteome_colors):
    this_cond_read_lengths = [this_read.query_length for val in cond_reads[this_proteome_cond].values()
                              for this_read in val]
    this_hist, bin_edges = np.histogram(this_cond_read_lengths, range=[0, 5000], bins=50)
    this_hist = this_hist / np.sum(this_hist)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    plt.plot(bin_centers, this_hist, c=this_color, label=f'{this_proteome_cond} [{len(this_cond_read_lengths)}]')
plt.legend(fontsize=8)
plt.xlabel('Read length (bps)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.subplot(num_rows, num_cols, 2)
for this_proteome_cond, this_color in zip(proteome_conditions, proteome_colors):
    this_cond_read_ids = [this_read.query_name for this_val in cond_reads[this_proteome_cond].values() for this_read in this_val]
    this_cond_polyA_lengths = collect_polyA_lengths_from_read_ids(df_polyA, this_cond_read_ids)
    this_hist, bin_edges = np.histogram(this_cond_polyA_lengths, range=[0, 500], bins=50)
    this_hist = this_hist / np.sum(this_hist)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    plt.plot(bin_centers, this_hist, c=this_color, label=f'{this_proteome_cond} [{len(this_cond_polyA_lengths)}]')
plt.legend(fontsize=8)
plt.xlabel('polyA length (bps)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)

for row_ind, num_bases in enumerate([5, 10, 20]):
    for col_ind, this_mod in enumerate(mods):
        for this_proteome_cond, this_color in zip(proteome_conditions, proteome_colors):
            plt.subplot(num_rows, num_cols, (row_ind+1)*num_cols+col_ind+1)
            all_gene_avg_mod_level = calc_mod_level_per_read_in_all_genes(cond_reads[this_proteome_cond], num_bases)
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
        plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $ per read, top {num_bases}', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)

plt.suptitle(f'{condition} {day}', fontsize=15)
plt.savefig(os.path.join(img_out, f'hist_avg_mods_per_read_polyA_len_grouped_by_protein_expression_{condition}_{day}.png'),
            bbox_inches='tight')
