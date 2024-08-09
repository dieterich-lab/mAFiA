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

conditions = ['SHAM', 'TAC']
day = 'day7'
num_bases = 5

proteome_conditions = [
    'up',
    'down'
]
proteome_colors = {
    'up': 'green',
    'down': 'red'
}

img_out = '/home/adrian/img_out/protein_abundance'
os.makedirs(img_out, exist_ok=True)

proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/protein_expression_data_mapped.tsv'
df_proteome = pd.read_csv(proteome_file, sep='\t')
sel_cols = list(df_proteome.columns[df_proteome.columns.str.contains(f'{day}')])
df_proteome = df_proteome[['gene'] + sel_cols]
for this_cond in conditions:
    this_cond_cols = list(df_proteome.columns[df_proteome.columns.str.contains(f'{this_cond}')])
    df_proteome[f'{this_cond}_{day}_avg'] = df_proteome[this_cond_cols].mean(axis=1)

thresh_percent_change = 0.05

df_proteome[f'percentage_change_{day}'] = (df_proteome[f'TAC_{day}_avg'] - df_proteome[f'SHAM_{day}_avg']) / df_proteome[f'SHAM_{day}_avg']
df_proteome_up = df_proteome[df_proteome[f'percentage_change_{day}'] >= thresh_percent_change].sort_values(f'SHAM_{day}_avg')
df_proteome_down = df_proteome[df_proteome[f'percentage_change_{day}'] < -thresh_percent_change].sort_values(f'SHAM_{day}_avg')

df_proteome_up.to_csv(os.path.join(img_out, f'protein_up_{day}.tsv'), sep='\t', index=False)
df_proteome_down.to_csv(os.path.join(img_out, f'protein_down_{day}.tsv'), sep='\t', index=False)

xlim = [18, 32]
ylim = [18, 32]

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(df_proteome[f'percentage_change_{day}'], bins=50, range=[-0.25, 0.25], log=True, fc='gray', alpha=0.5)
plt.hist(df_proteome_up[f'percentage_change_{day}'], bins=50, range=[-0.25, 0.25], log=True, fc=proteome_colors['up'])
plt.hist(df_proteome_down[f'percentage_change_{day}'], bins=50, range=[-0.25, 0.25], log=True, fc=proteome_colors['down'])
plt.xlabel('% change', fontsize=12)
plt.ylabel('Protein count', fontsize=12)
plt.subplot(1, 2, 2)
plt.scatter(df_proteome[f'SHAM_{day}_avg'], df_proteome[f'TAC_{day}_avg'], s=0.5, c='gray', alpha=0.5)
plt.scatter(df_proteome_up[f'SHAM_{day}_avg'], df_proteome_up[f'TAC_{day}_avg'], s=5, c=proteome_colors['up'])
plt.scatter(df_proteome_down[f'SHAM_{day}_avg'], df_proteome_down[f'TAC_{day}_avg'], s=5, c=proteome_colors['down'])
num_up = len(df_proteome_up)
x_pos_up = np.linspace(xlim[0], xlim[1]*0.9, num_up)
# y_pos_up = np.linspace(ylim[0], ylim[1], num_up)
# y_pos_up = y_pos_up + 0.1
y_pos_up = ylim[1] - np.arange(num_up)[::-1] * 0.5
num_down = len(df_proteome_down)
x_pos_down = np.linspace(xlim[0]*1.1, xlim[1], num_down)
# y_pos_down = np.linspace(ylim[0], ylim[1], num_down)
# y_pos_down = y_pos_down - 0.1
y_pos_down = ylim[0] + np.arange(num_down) * 0.5
for ind, (_, this_row) in enumerate(df_proteome_up.iterrows()):
    (x, y) = (this_row[f'SHAM_{day}_avg'], this_row[f'TAC_{day}_avg'])
    (x_text, y_text) = (x_pos_up[ind], y_pos_up[ind])
    plt.annotate(this_row['gene'], (x, y), (x_text, y_text), fontsize=8,
                 c=proteome_colors['up'], arrowprops=dict(facecolor=proteome_colors['up'], width=0.1, headwidth=0.1))
for ind, (_, this_row) in enumerate(df_proteome_down.iterrows()):
    (x, y) = (this_row[f'SHAM_{day}_avg'], this_row[f'TAC_{day}_avg'])
    (x_text, y_text) = (x_pos_down[ind], y_pos_down[ind])
    plt.annotate(this_row['gene'], (x, y), (x_text, y_text), fontsize=8,
                 c=proteome_colors['down'], arrowprops=dict(facecolor=proteome_colors['down'], width=0.1, headwidth=0.1))
plt.xlabel('SHAM', fontsize=12)
plt.ylabel('TAC', fontsize=12)
plt.title('Protein Abundance', fontsize=12)
plt.suptitle(day, fontsize=15)
plt.savefig(os.path.join(img_out, f'protein_abundance_{day}.png'), bbox_inches='tight')

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

bam_files = {
    this_cond: f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/{this_cond}_{day}/chrALL.mAFiA.reads.bam'
    for this_cond in conditions
}
# polyA_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA/{condition}_{day}_read2polyA_length.txt'

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


def calc_mod_level_per_read(in_reads, in_num_bases):
    read_avg_mod_level = {}
    for this_read in in_reads:
        read_avg_mod_level[this_read.query_name] = {}
        for mod in dict_mod_code.keys():
            mod_bases = this_read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])
            if len(mod_bases) >= in_num_bases:
                read_avg_mod_level[this_read.query_name][mod] = np.mean(
                    np.partition([tup[1] for tup in mod_bases], -in_num_bases)[-in_num_bases:]) / 255.0
            else:
                read_avg_mod_level[this_read.query_name][mod] = np.mean([tup[1] for tup in mod_bases]) / 255.0
    return read_avg_mod_level


def calc_mod_level_per_read_in_all_genes(in_gene_reads, in_num_bases):
    gene_avg_mod_level = {}
    for this_gene, this_gene_reads in tqdm(in_gene_reads.items()):
        this_gene_read_avg_mod_level = calc_mod_level_per_read(this_gene_reads, in_num_bases)
        if len(this_gene_read_avg_mod_level):
            gene_avg_mod_level[this_gene] = this_gene_read_avg_mod_level
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

thresh_percent_change_general = 0.02
df_proteome_up_general = df_proteome[df_proteome[f'percentage_change_{day}'] >= thresh_percent_change_general].sort_values(f'SHAM_{day}_avg')
df_proteome_down_general = df_proteome[df_proteome[f'percentage_change_{day}'] < -thresh_percent_change_general].sort_values(f'SHAM_{day}_avg')

cond_df_proteome = {
    'up': df_proteome_up_general,
    'down': df_proteome_down_general
}

cond_exp_reads = {
    this_cond: {
        this_exp: collect_proteome_reads(this_df, bam_files[this_exp])
        for this_exp in ['SHAM', 'TAC']
    }
    for this_cond, this_df in cond_df_proteome.items()
}

num_bins = 10
max_bin = 1.0

num_rows = 2
num_cols = 2
num_bases = 5
experiments = ['SHAM', 'TAC']
exp_color = {
    'SHAM': 'skyblue',
    'TAC': 'midnightblue'
}

hist_range = [0, 1]
num_bins = 20
bin_edges = np.linspace(*hist_range, num_bins+1)
bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
ylim = [0, 0.12]

plt.figure(figsize=(10, 10))
for row_ind, this_exp in enumerate(experiments):
    for this_proteome_cond in ['up', 'down']:
        exp_gene_avg_mod = calc_mod_level_per_read_in_all_genes(cond_exp_reads[this_proteome_cond][this_exp], num_bases)
        for col_ind, this_mod in enumerate(mods):
            this_mod_vals = [this_read_mods.get(this_mod, np.nan)
                             for val in exp_gene_avg_mod.values()
                             for this_read_mods in val.values()]
            plt.subplot(2, 2, row_ind*2+col_ind+1)
            this_hist, _ = np.histogram(this_mod_vals, range=hist_range, bins=num_bins)
            this_hist = this_hist / np.sum(this_hist)
            plt.plot(bin_centers, this_hist, c=proteome_colors[this_proteome_cond], label=this_proteome_cond)
            # plt.hist(this_mod_vals, range=[0, 1], bins=20,
            #          fc=proteome_colors[this_proteome_cond], label=this_proteome_cond, alpha=0.5)
            plt.legend()
            plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $ per read, top {num_bases}',
                       fontsize=12)
            plt.ylabel('Frequency', fontsize=12)
            # plt.ylim(ylim)
            plt.title(this_exp, fontsize=12)
plt.suptitle(day, fontsize=15)
plt.savefig(os.path.join(img_out, f'avg_mod_per_read_split_by_protein_changes_{day}.png'), bbox_inches='tight')

# hist_range = [-0.5, 0.5]
# num_bins = 10
# bin_edges = np.linspace(*hist_range, num_bins+1)
# bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
#
# plt.figure(figsize=(10, 4))
# for this_proteome_cond in ['down', 'up']:
#     gene_mod_change = {}
#     for this_gene in cond_df_proteome[this_proteome_cond]['gene']:
#         if (this_gene in cond_exp_reads[this_proteome_cond]['SHAM'].keys()) and (this_gene in cond_exp_reads[this_proteome_cond]['SHAM'].keys()):
#             gene_mod_change[this_gene] = {}
#             tac_read_avg_mod = calc_mod_level_per_read(cond_exp_reads[this_proteome_cond]['TAC'][this_gene], num_bases)
#             sham_read_avg_mod = calc_mod_level_per_read(cond_exp_reads[this_proteome_cond]['SHAM'][this_gene], num_bases)
#             tac_m6A_mean = np.nanmean([this_val.get('m6A', np.nan) for this_val in tac_read_avg_mod.values()])
#             tac_psi_mean = np.nanmean([this_val.get('psi', np.nan) for this_val in tac_read_avg_mod.values()])
#             sham_m6A_mean = np.nanmean([this_val.get('m6A', np.nan) for this_val in sham_read_avg_mod.values()])
#             sham_psi_mean = np.nanmean([this_val.get('psi', np.nan) for this_val in sham_read_avg_mod.values()])
#             gene_mod_change[this_gene]['m6A'] = (tac_m6A_mean - sham_m6A_mean) / sham_m6A_mean
#             gene_mod_change[this_gene]['psi'] = (tac_psi_mean - sham_psi_mean) / sham_psi_mean
#     for mod_ind, this_mod in enumerate(mods):
#         this_mod_keys = list(gene_mod_change.keys())
#         this_mod_vals = [v[this_mod] for v in gene_mod_change.values()]
#         protein_deltas = [cond_df_proteome[this_proteome_cond][cond_df_proteome[this_proteome_cond]['gene']==this_gene][f'percentage_change_{day}'].values[0]
#                           for this_gene in this_mod_keys]
#         # this_hist, _ = np.histogram(this_mod_vals, range=hist_range, bins=num_bins)
#         # this_hist = this_hist / np.sum(this_hist)
#
#         # mask = ~np.isnan(this_mod_vals)
#
#         plt.subplot(1, 2, mod_ind+1)
#         plt.scatter(this_mod_vals, protein_deltas, c=proteome_colors[this_proteome_cond], s=2, label=this_proteome_cond)
#         plt.legend()
#         plt.xlabel(rf'$\Delta \langle P({{{dict_mod_display[this_mod]}}}) \rangle$',
#                    fontsize=12)
#         plt.ylabel('% change protein', fontsize=12)