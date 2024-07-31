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

dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

day = 'day7'
polyA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA'
polyA_sham_file = os.path.join(polyA_dir, f'SHAM_{day}_read2polyA_length.txt')
polyA_tac_file = os.path.join(polyA_dir, f'TAC_{day}_read2polyA_length.txt')
sham_bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_{day}/chrALL.mAFiA.reads.bam'
tac_bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_{day}/chrALL.mAFiA.reads.bam'


img_out = '/home/adrian/img_out/polyA_len'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

df_polyA_sham = pd.read_csv(polyA_sham_file, sep='\t', names=['read_id', 'polyA_len'])
df_polyA_tac = pd.read_csv(polyA_tac_file, sep='\t', names=['read_id', 'polyA_len'])

dict_polyA_sham = {k: v for k, v in df_polyA_sham[['read_id', 'polyA_len']].values}
dict_polyA_tac = {k: v for k, v in df_polyA_tac[['read_id', 'polyA_len']].values}

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.hist(dict_polyA_sham.values(), range=[0, 400], bins=100)
plt.subplot(1, 2, 2)
plt.hist(dict_polyA_tac.values(), range=[0, 400], bins=100)

# thresh_polyA_len = 90
upper_lim = 200
lower_lim = 50
sham_above_read_ids = df_polyA_sham[df_polyA_sham['polyA_len'] >= upper_lim]['read_id'].tolist()
sham_below_read_ids = df_polyA_sham[df_polyA_sham['polyA_len'] < lower_lim]['read_id'].tolist()
with open(os.path.join(polyA_dir, f'SHAM_{day}_read_ids_polyA_len_above_{upper_lim}.txt'), 'w') as f_out:
    for this_read_id in sham_above_read_ids:
        f_out.writelines(this_read_id + '\n')
# with open(os.path.join(polyA_dir, f'SHAM_{day}_read_ids_polyA_len_below_{lower_lim}.txt'), 'w') as f_out:
#     for this_read_id in sham_below_read_ids:
#         f_out.writelines(this_read_id + '\n')

# sham_read_ids = list(df_polyA_sham['read_id'].values)
# sham_read_percentage_mods = {}
#
# with pysam.AlignmentFile(sham_bam_file, 'rb') as sham_bam:
#     for this_read in tqdm(sham_bam.fetch()):
#         if this_read.query_name in sham_read_ids:
#             dict_base_counts = Counter(this_read.query_sequence)
#             sham_read_percentage_mods[this_read.query_name] = {}
#             for this_mod in dict_mod_code.keys():
#                 mod_bases = this_read.modified_bases.get(('N', 0, dict_mod_code[this_mod]), [])
#                 if this_mod=='m6A':
#                     base_count = dict_base_counts['A']
#                 elif this_mod=='psi':
#                     base_count = dict_base_counts['T']
#                 sham_read_percentage_mods[this_read.query_name][this_mod] = np.sum([this_mod_base[1]>=128
#                                                                              for this_mod_base in mod_bases]) / base_count
#
# mods_tail_len = [(v['m6A'], v['psi'], dict_polyA_sham[k]) for k, v in sham_read_percentage_mods.items()]
# vec_m6A, vec_psi, vec_tail_len = np.vstack(mods_tail_len).T
#
# plt.figure(figsize=(10, 5))
# plt.subplot(1, 2, 1)
# plt.scatter(vec_m6A, vec_tail_len, s=0.5)
# plt.xlabel(rf"$\%{{{dict_mod_display['m6A']}}}$ in read", fontsize=12)
# plt.ylabel('polyA tail length (bps)', fontsize=12)
# plt.subplot(1, 2, 2)
# plt.scatter(vec_psi, vec_tail_len, s=0.5)
# plt.xlabel(rf"$\%{{{dict_mod_display['psi']}}}$ in read", fontsize=12)
# plt.ylabel('polyA tail length (bps)', fontsize=12)
#
# num_bins = 10
# bin_width = 0.02
# bin_edges = np.arange(0.0, (num_bins+1)*bin_width, bin_width)
# mat_tail_len = np.zeros((num_bins, num_bins))
# for i_ind in range(len(bin_edges)-1):
#     i_start = bin_edges[i_ind]
#     i_end = bin_edges[i_ind+1]
#     i_mask = (vec_psi>=i_start) * (vec_psi<i_end)
#     for j_ind in range(len(bin_edges)-1):
#         j_start = bin_edges[j_ind]
#         j_end = bin_edges[j_ind + 1]
#         j_mask = (vec_m6A>=j_start) * (vec_m6A<j_end)
#         mat_tail_len[i_ind, j_ind] = np.mean(vec_tail_len[i_mask * j_mask])
#
# plt.figure(figsize=(5, 5))
# im = plt.imshow(mat_tail_len, origin='lower')
# plt.xticks(np.arange(num_bins+1)-0.5, bin_edges)
# plt.yticks(np.arange(num_bins+1)-0.5, bin_edges)
# plt.xlabel(rf"$\%{{{dict_mod_display['m6A']}}}$ in read", fontsize=12)
# plt.ylabel(rf"$\%{{{dict_mod_display['psi']}}}$ in read", fontsize=12)
# clb = plt.colorbar(im)