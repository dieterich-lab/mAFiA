import os
import pandas as pd
import numpy as np
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def get_single_read_mod_bases(in_bam_file):
    read_mod_bases = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for this_read in tqdm(bam.fetch()):
            dict_read_ref = {tup[0]: tup[1] for tup in this_read.get_aligned_pairs(matches_only=True)}
            # int_strand = 0 if this_read.flag==0 else 1

            this_read_mod_bases = {
                k: {} for k, v in mod_code.items()
            }
            for k, v in mod_code.items():
                for (read_pos, score) in this_read.modified_bases.get(('N', 0, v), []):
                    this_read_mod_bases[k][dict_read_ref[read_pos]] = score

            read_mod_bases.append(this_read_mod_bases)
    return read_mod_bases

def get_min_max_pos_in_read(in_read_mod_bases):
    out_read_mod_bases = []
    for this_read_mod_bases in tqdm(in_read_mod_bases):
        all_positions = [pos for this_list in this_read_mod_bases.values() for pos in list(this_list.keys())]
        if len(all_positions):
            this_read_mod_bases['min'] = np.min(all_positions)
            this_read_mod_bases['max'] = np.max(all_positions)
            out_read_mod_bases.append(this_read_mod_bases)
    return out_read_mod_bases

def get_position_matched_scores(read_mod_bases, m6A_pos, psi_pos):
    candidate_mod_bases = [this_mod_base for this_mod_base in read_mod_bases if
                           min(m6A_pos, psi_pos)>=this_mod_base['min']
                           and max(m6A_pos, psi_pos)<this_mod_base['max']
                           ]
    pos_matched_scores = []
    for this_read_mod_bases in candidate_mod_bases:
        if (m6A_pos in this_read_mod_bases['m6A'].keys()) and (psi_pos in this_read_mod_bases['psi'].keys()):
            pos_matched_scores.append([this_read_mod_bases['m6A'][m6A_pos], this_read_mod_bases['psi'][psi_pos]])
    return pos_matched_scores

chrom = 'X'
workspace = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT/chr{chrom}'
bed_file = os.path.join(workspace, 'mAFiA.sites.bed')
bam_file = os.path.join(workspace, 'mAFiA.reads.bam')

img_out = '/home/adrian/img_out/psi-mAFiA_cross_talk'
os.makedirs(img_out, exist_ok=True)

thresh_conf = 80.0
thresh_modRatio = 0.0

df_sites = pd.read_csv(bed_file, sep='\t', dtype={'chrom': str})
df_sites = df_sites[
    (df_sites['confidence']>=thresh_conf)
    * (df_sites['modRatio']>=thresh_modRatio)
    ]
psi_positions = df_sites[df_sites['name']=='psi']['chromStart'].values
m6A_positions = df_sites[df_sites['name']=='m6A']['chromStart'].values

mod_code = {
    'm6A': 21891,
    'psi': 17802
}

### read-level correlation ###
single_read_mod_bases = get_single_read_mod_bases(bam_file)
filtered_read_mod_bases = get_min_max_pos_in_read(single_read_mod_bases)

pos_pairs = []
vec_corr = []
vec_r = []
for this_psi_pos in tqdm(psi_positions):
    for this_m6A_pos in m6A_positions:
        dist = abs(this_psi_pos-this_m6A_pos)
        if dist>=10000:
            continue
        read_pos_matched_scores = get_position_matched_scores(filtered_read_mod_bases, this_m6A_pos, this_psi_pos)
        if len(read_pos_matched_scores):
            corr = np.corrcoef(np.vstack(read_pos_matched_scores).T)[0, 1]
            if not np.isnan(corr):
                pos_pairs.append((this_psi_pos, this_m6A_pos))
                vec_corr.append(corr)
                vec_r.append(dist)

vec_corr = np.array(vec_corr)
vec_r = np.array(vec_r)
pos_psi = [tup[0] for tup in pos_pairs]
pos_m6A = [tup[1] for tup in pos_pairs]
df_read_level = pd.DataFrame({'chromStart_psi': pos_psi, 'chromStart_m6A': pos_m6A, 'dist': vec_r, 'corr': vec_corr})

### plot ###
x_min = 0
x_max = 10000
num_bins = 20
# bin_x = np.logspace(1, 4, num_bins+1)
bin_x = np.linspace(x_min, x_max, num_bins+1)
xticks = np.linspace(x_min, x_max, 6)

mean_corr = []
mean_abs_corr = []
for i_start in range(num_bins):
    r_start = bin_x[i_start]
    r_end = bin_x[i_start+1]
    this_bin_corr = vec_corr[(vec_r>=r_start) * (vec_r<r_end)]
    mean_corr.append(np.nanmean(this_bin_corr))
    mean_abs_corr.append(np.nanmean(np.abs(this_bin_corr)))

bin_center_x = 0.5 * (bin_x[1:] + bin_x[:-1])

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(vec_r, vec_corr, '.')
# plt.xscale('log')
plt.yticks(np.linspace(-1, 1, 5))
plt.xlabel('Distance (nts)', fontsize=12)
plt.ylabel('Read-level correlation C(r)', fontsize=12)
plt.subplot(1, 2, 2)
plt.plot(bin_center_x, mean_corr, 'b-', label=r'$\langle C(r) \rangle$')
plt.plot(bin_center_x, mean_abs_corr, 'r-', label=r'$\langle\left| C(r) \right|\rangle$')
# plt.yticks(np.linspace(-1, 1, 5))
# plt.xticks([100, 1000, 10000])
plt.yticks(np.arange(0, 0.21, 0.1))
plt.legend(loc='upper right', fontsize=10)
plt.xlabel('Distance (nts)', fontsize=12)
plt.ylabel(r'Avg. read-level correlation $\langle C(r) \rangle$', fontsize=12)
plt.suptitle(f'HEK293 chr{chrom}', fontsize=15)
plt.savefig(os.path.join(img_out, f'read_level_correlation_chr{chrom}.png'), bbox_inches='tight')

### site-level correlation ###
site_level_corr = []
for start_i in tqdm(range(len(bin_x)-1)):
    start_dist = bin_x[start_i]
    end_dist = bin_x[start_i+1]
    sub_df = df_read_level[(df_read_level['dist']>=start_dist) * (df_read_level['dist']<end_dist)]
    vec_psi_modRatio = [df_sites[(df_sites['chromStart']==pos) * (df_sites['name']=='psi')]['modRatio'].values[0] for pos in sub_df['chromStart_psi'].values]
    vec_m6A_modRatio = [df_sites[(df_sites['chromStart']==pos) * (df_sites['name']=='m6A')]['modRatio'].values[0] for pos in sub_df['chromStart_m6A'].values]
    site_level_corr.append(np.corrcoef(np.vstack([vec_psi_modRatio, vec_m6A_modRatio]))[0, 1])

plt.figure(figsize=(5, 5))
plt.plot(bin_center_x, site_level_corr)
plt.yticks(np.linspace(-0.2, 0.2, 5))
plt.xlabel('Distance (nts)', fontsize=12)
plt.ylabel('Site-level correlation', fontsize=12)
plt.suptitle(f'HEK293 chr{chrom}', fontsize=15)
plt.savefig(os.path.join(img_out, f'site_level_correlation_chr{chrom}.png'), bbox_inches='tight')
