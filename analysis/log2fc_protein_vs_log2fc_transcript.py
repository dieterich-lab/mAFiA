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

proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'
df_proteome = pd.read_csv(proteome_file, sep='\t')
df_proteome_thresh = df_proteome[
    (df_proteome['tp'] == 'Combined')
    # * (df_proteome['p.adj'] < thresh_pval)
]
df_proteome_thresh.rename(columns={'logFC': 'log2fc_proteome'}, inplace=True)

transcriptome_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_abundance/TAC_TAC_merged_vs_SHAM_merged.tsv'
df_transcriptome = pd.read_csv(transcriptome_file, sep='\t')
df_transcriptome.rename(columns={'log2fc': 'log2fc_transcriptome'}, inplace=True)

gene_transcriptome_proteome = {}
for this_gene in df_proteome_thresh['SYMBOL']:
    if this_gene in df_transcriptome['gene'].values:
        gene_transcriptome_proteome[this_gene] = (
            df_transcriptome[df_transcriptome['gene'] == this_gene]['log2fc_transcriptome'].values[0],
            df_proteome_thresh[df_proteome_thresh['SYMBOL'] == this_gene]['log2fc_proteome'].values[0]
        )

vec_x, vec_y = np.vstack(gene_transcriptome_proteome.values()).T

plt.figure(figsize=(5, 5))
plt.scatter(vec_x, vec_y, s=1)
plt.xlim([-3, 3])
plt.ylim([-1, 1])
plt.xlabel('log2fc transcript')
plt.ylabel('log2fc protein')