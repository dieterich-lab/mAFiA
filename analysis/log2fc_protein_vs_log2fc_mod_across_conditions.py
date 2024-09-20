import os
import pandas as pd
import re
import numpy as np
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}


def get_gene_avg_mod_prob(in_bam, in_row, thresh_coverage=10):
    flag_required = 0 if in_row['strand'] == '+' else 16
    mod_read_mod_probs = {this_mod: [] for this_mod in mods}
    read_count = 0
    chrom, chromStart, chromEnd = in_row[['chrom', 'chromStart', 'chromEnd']]
    for this_read in in_bam.fetch(chrom, chromStart, chromEnd):
        if this_read.flag == flag_required:
            query_to_ref_pos = {pair[0]: pair[1] for pair in this_read.get_aligned_pairs(matches_only=True)}
            for mod in mods:
                ref_pos_mod_prob = [(query_to_ref_pos[pair[0]], pair[1])
                                    for pair in this_read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])
                                    ]
                if len(ref_pos_mod_prob):
                    read_count += 1
                    mod_read_mod_probs[mod].extend([
                        pair for pair in ref_pos_mod_prob
                        if (pair[0] >= chromStart) and (pair[0] < chromEnd)
                    ])
    # print(read_count)
    if read_count >= thresh_coverage:
        return {mod: np.mean([pair[1] >= 128 for pair in mod_read_mod_probs[mod]]) for mod in mods}
    else:
        return {mod: np.nan for mod in mods}


thresh_pval = 0.05

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
cond_names = [this_cond.rstrip('_merged') for this_cond in conditions]

img_out = '/home/adrian/img_out/protein_expression'
os.makedirs(img_out, exist_ok=True)

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

protein_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'
df_protein = pd.read_csv(protein_file, sep='\t')
df_protein_thresh = df_protein[
    (df_protein['tp'] == 'Combined')
    * (df_protein['p.adj'] < thresh_pval)
]
df_protein_thresh.rename(columns={'logFC': 'log2fc_protein'}, inplace=True)

gene_bed_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
gene_bed_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
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
df_gene = pd.read_csv(gene_bed_file, sep='\t', names=bed_fields)
df_gene['gene_name'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

gene_cond_avg_mod_prob = {}
for this_gene in tqdm(df_protein_thresh['SYMBOL']):
    # print(this_gene)
    sub_df_gene = df_gene[df_gene['gene_name'] == this_gene]
    if len(sub_df_gene) == 1:
        gene_cond_avg_mod_prob[this_gene] = {}
        for this_cond in conditions:
            this_cond_bam_path = os.path.join(res_dir, ds, this_cond, 'chrALL.mAFiA.reads.bam')
            with pysam.AlignmentFile(this_cond_bam_path, 'rb') as this_cond_bam:
                gene_cond_avg_mod_prob[this_gene][this_cond] = get_gene_avg_mod_prob(this_cond_bam, sub_df_gene.iloc[0])
    elif len(sub_df_gene) > 1:
        print('Duplicate annotation!', sub_df_gene)

gene_log2fc_mod = {}
for this_gene in df_protein_thresh['SYMBOL']:
    gene_log2fc_mod[this_gene] = {}
    if this_gene in gene_cond_avg_mod_prob.keys():
        for this_mod in mods:
            gene_log2fc_mod[this_gene][this_mod] = np.log2(gene_cond_avg_mod_prob[this_gene][conditions[1]][this_mod]
                                                           / gene_cond_avg_mod_prob[this_gene][conditions[0]][this_mod])
    else:
        gene_log2fc_mod[this_gene] = {this_mod: np.nan for this_mod in mods}

for this_mod in mods:
    df_protein_thresh[f'log2fc_{this_mod}'] = [gene_log2fc_mod[this_gene][this_mod] for this_gene in df_protein_thresh['SYMBOL']]

df_protein_thresh.sort_values('log2fc_protein', inplace=True)
print(df_protein_thresh[['log2fc_protein', 'log2fc_m6A', 'log2fc_psi']])

df_protein_thresh.to_csv(os.path.join(img_out, 'df_protein_expression_with_log2fc_mods.tsv'), sep='\t', index=False)

boundary = 0.05

xmax = 1.0
ymax = 0.06
xticks = np.round(np.linspace(-xmax, xmax, 5), 2)
yticks = np.round(np.linspace(-ymax, ymax, 5), 2)

num_bins = 4
bin_edges = np.round(np.linspace(-xmax, xmax, num_bins+1), 2)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    binned_y = []
    vec_x = df_protein_thresh[f'log2fc_protein'].values
    vec_y = df_protein_thresh[f'log2fc_{this_mod}'].values
    # binned_y = [
    #      [x for x in df_protein_thresh[f'log2fc_protein'].values[(vec_x < -boundary) * (vec_x >= -1)] if ~np.isnan(x) and ~np.isinf(x)],
    #     [x for x in df_protein_thresh[f'log2fc_protein'].values[(vec_x >= boundary) * (vec_x < 1)] if ~np.isnan(x) and ~np.isinf(x)],
    # ]
    for this_bin_ind in range(len(bin_edges)-1):
        bin_start = bin_edges[this_bin_ind]
        bin_end = bin_edges[this_bin_ind+1]
        mask = (vec_x >= bin_start) * (vec_x < bin_end)
        binned_y.append([x for x in vec_y[mask] if ~np.isnan(x)])
    print(this_mod, [len(ll) for ll in binned_y])
    # plt.boxplot(binned_y, positions=[-0.5, 0.5], widths=0.25, whis=0.5, showfliers=False)
    # plt.boxplot(binned_y, positions=bin_centers, widths=xmax/len(bin_centers), whis=1.5, showfliers=False)

    median_binned_y = [np.median(this_bin) for this_bin in binned_y]
    plt.bar(bin_centers, median_binned_y, width=xmax/len(bin_centers))

    plt.xlim([-xmax, xmax])
    plt.ylim([0, ymax*1.05])
    # plt.xticks([-0.5, 0.5], [f'< -{boundary} ({len(binned_y[0])})', rf'$\geq$ {boundary} ({len(binned_y[1])})'])
    # plt.xticks(xticks)
    plt.xticks(bin_edges, bin_edges)
    plt.yticks(np.linspace(0, ymax, 5))
    # plt.yticks(yticks)
    plt.xlabel(r'$log_{2}fc$ protein')
    plt.ylabel(r'$\langle$ $log_{2}fc$ $\bar{S}$ $\rangle$')
    plt.title(rf'${{{dict_mod_display[this_mod]}}}$')
plt.suptitle(f'{conditions[1]} vs {conditions[0]}')
plt.savefig(os.path.join(img_out, 'log2fc_mod_vs_log2fc_protein.png'), bbox_inches='tight')
