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


thresh_pval = 0.05

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
cond_names = [this_cond.rstrip('_merged') for this_cond in conditions]

img_out = '/home/adrian/img_out/protein_expression'
os.makedirs(img_out, exist_ok=True)

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'
df_proteome_all = pd.read_csv(proteome_file, sep='\t')
df_proteome = df_proteome_all[
    (df_proteome_all['tp'] == 'Combined')
    * (df_proteome_all['p.adj'] < thresh_pval)
]
df_proteome.rename(columns={'logFC': 'log2fc_proteome'}, inplace=True)

epitranscriptome_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_stoichiometry/TAC_TAC_merged_vs_SHAM_merged.tsv'
df_epitranscriptome_all = pd.read_csv(epitranscriptome_file, sep='\t')
df_epitranscriptome = df_epitranscriptome_all[
    df_epitranscriptome_all['pval'] < thresh_pval
]
# df_epitranscriptome = df_epitranscriptome_all
# df_epitranscriptome = df_epitranscriptome_all[
#     (df_epitranscriptome_all['num_nts_SHAM_merged'] >= 1000)
#     * (df_epitranscriptome_all['num_nts_TAC_merged'] >= 1000)
# ]
df_epitranscriptome.rename(columns={'log2fc': 'log2fc_epitranscriptome'}, inplace=True)

transcriptome_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_abundance/TAC_TAC_merged_vs_SHAM_merged.tsv'
df_transcriptome = pd.read_csv(transcriptome_file, sep='\t')
df_transcriptome.rename(columns={'log2fc': 'log2fc_transcriptome'}, inplace=True)

gene_transcriptome_epitranscriptome_proteome = {}
for this_gene in tqdm(df_proteome['SYMBOL']):
    # if (this_gene in df_transcriptome['gene'].values) and (this_gene in df_epitranscriptome['gene'].values):
    if (this_gene in df_transcriptome['gene'].values):
        sub_df_epitranscriptome = df_epitranscriptome[df_epitranscriptome['gene'] == this_gene]
        log2fc_mods = {}
        for this_mod in mods:
            if this_mod in sub_df_epitranscriptome['mod'].values:
                log2fc_mods[this_mod] = sub_df_epitranscriptome[
                    sub_df_epitranscriptome['mod'] == this_mod
                    ]['log2fc_epitranscriptome'].values[0]
            else:
                log2fc_mods[this_mod] = None
        gene_transcriptome_epitranscriptome_proteome[this_gene] = (
            df_transcriptome[df_transcriptome['gene'] == this_gene]['log2fc_transcriptome'].values[0],
            log2fc_mods['m6A'],
            log2fc_mods['psi'],
            df_proteome[df_proteome['SYMBOL'] == this_gene]['log2fc_proteome'].values[0]
        )

vec_abundance, vec_m6A, vec_psi, vec_protein = np.vstack(gene_transcriptome_epitranscriptome_proteome.values()).T

### scatter plots ###
# plt.figure(figsize=(10, 10))
# plt.subplot(2, 2, 1)
# plt.scatter(vec_abundance, vec_protein, s=1)
# # plt.xlim([-3, 3])
# # plt.ylim([-1, 1])
# plt.xlabel('log2fc transcript abundance')
# plt.ylabel('log2fc protein')
# plt.subplot(2, 2, 2)
# plt.scatter(vec_m6A, vec_psi, s=1)
# # plt.xlim([-3, 3])
# # plt.ylim([-1, 1])
# plt.xlabel('log2fc m6A')
# plt.ylabel('log2fc psi')
# plt.subplot(2, 2, 3)
# plt.scatter(vec_m6A, vec_protein, s=1)
# # plt.xlim([-3, 3])
# # plt.ylim([-1, 1])
# plt.xlabel('log2fc m6A')
# plt.ylabel('log2fc protein')
# plt.subplot(2, 2, 4)
# plt.scatter(vec_psi, vec_protein, s=1)
# # plt.xlim([-3, 3])
# # plt.ylim([-1, 1])
# plt.xlabel('log2fc psi')
# plt.ylabel('log2fc protein')

### boxplot ###
def get_1d_binned_values(in_vec_x, in_vec_y, xmax, num_bins):
    in_vec_x[in_vec_x == None] = -np.inf

    x_bin_edges = np.linspace(-xmax, xmax, num_bins + 1)
    x_bin_centers = 0.5 * (x_bin_edges[1:] + x_bin_edges[:-1])
    binned_y = []
    for bin_i in range(len(x_bin_edges) - 1):
        bin_start = x_bin_edges[bin_i]
        bin_stop = x_bin_edges[bin_i + 1]
        mask = (in_vec_x >= bin_start) * (in_vec_x < bin_stop)
        binned_y.append(in_vec_y[mask])
    medians_y = [np.mean(ll) for ll in binned_y]
    return x_bin_edges, x_bin_centers, binned_y, medians_y

ymax = 0.55
yticks = np.linspace(-0.5, 0.5, 5)

plt.figure(figsize=(15, 4))

plt.subplot(1, 3, 1)
xmax = 2
x_bin_edges, x_bin_centers, binned_y, medians_y = get_1d_binned_values(vec_abundance, vec_protein, xmax, 4)
plt.boxplot(binned_y, positions=x_bin_centers, showfliers=False, whis=0.5, medianprops={'linewidth': 0})
plt.plot(x_bin_centers, medians_y, 'r-o')
plt.xticks(x_bin_edges, x_bin_edges)
plt.yticks(yticks)
plt.xlim([-xmax, xmax])
plt.ylim([-ymax, ymax])
plt.xlabel('log2fc transcript abundance')
plt.ylabel('log2fc protein abundance')

plt.subplot(1, 3, 2)
xmax = 1
x_bin_edges, x_bin_centers, binned_y, medians_y = get_1d_binned_values(vec_m6A, vec_protein, xmax, 4)
plt.boxplot(binned_y, positions=x_bin_centers, showfliers=False, whis=0.5, medianprops={'linewidth': 0})
plt.plot(x_bin_centers, medians_y, 'r-o')
plt.xticks(x_bin_edges, x_bin_edges)
plt.yticks(yticks)
plt.xlim([-xmax, xmax])
plt.ylim([-ymax, ymax])
plt.xlabel(f"log2fc transcript ${{{dict_mod_display['m6A']}}}$")
plt.ylabel('log2fc protein abundance')

plt.subplot(1, 3, 3)
xmax = 1
x_bin_edges, x_bin_centers, binned_y, medians_y = get_1d_binned_values(vec_psi, vec_protein, xmax, 4)
plt.boxplot(binned_y, positions=x_bin_centers, showfliers=False, whis=0.5, medianprops={'linewidth': 0})
plt.plot(x_bin_centers, medians_y, 'r-o')
plt.xticks(x_bin_edges, x_bin_edges)
plt.yticks(yticks)
plt.xlim([-xmax, xmax])
plt.ylim([-ymax, ymax])
plt.xlabel(f"log2fc transcript ${{{dict_mod_display['psi']}}}$")
plt.ylabel('log2fc protein abundance')

plt.savefig(os.path.join(img_out, 'log2fc_protein_vs_transcript_abundance_stoichiometry.png'), bbox_inches='tight')


### 2d heatmap ###
# xmax = 0.5
# ymax = 0.5
# # bin_width = 0.5
# num_bins = 2
# x_bin_edges = np.linspace(-xmax, xmax, num_bins+1)
# y_bin_edges = np.linspace(-ymax, ymax, num_bins+1)
#
# vec_mods = {
#     'm6A': vec_m6A,
#     'psi': vec_psi
# }
#
# mod_binned_z = {}
# for this_mod in mods:
#     binned_z = [[[] for j in range(len(x_bin_edges)-1)] for i in range(len(y_bin_edges)-1)]
#     for bin_i in range(len(y_bin_edges)-1):
#         y_bin_start = y_bin_edges[bin_i]
#         y_bin_stop = y_bin_edges[bin_i+1]
#         y_mask = (vec_mods[this_mod] >= y_bin_start) * (vec_mods[this_mod] < y_bin_stop)
#         for bin_j in range(len(x_bin_edges)-1):
#             x_bin_start = x_bin_edges[bin_j]
#             x_bin_stop = x_bin_edges[bin_j+1]
#             x_mask = (vec_abundance >= x_bin_start) * (vec_abundance < x_bin_stop)
#             binned_z[bin_i][bin_j] = np.median(vec_protein[x_mask * y_mask])
#     mod_binned_z[this_mod] = np.array(binned_z)
#
# plt.figure(figsize=(10, 5))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     plt.imshow(mod_binned_z[this_mod], origin='lower')
#     plt.xticks(np.arange(len(x_bin_edges))-0.5, x_bin_edges)
#     plt.yticks(np.arange(len(y_bin_edges))-0.5, y_bin_edges)
#     plt.xlabel('log2fc transcript abundance')
#     plt.ylabel('log2fc stoichiometry')
#     plt.title(f'${{{dict_mod_display[this_mod]}}}$')