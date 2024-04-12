import os
import pandas as pd
import numpy as np
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

mod_code = {
    'm6A': 21891,
    'psi': 17802
}

def get_read_mod_bases(in_bam_file):
    read_mod_bases = {}
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for this_read in tqdm(bam.fetch()):
            this_read_mod_bases = {}
            for k, v in mod_code.items():
                this_read_mod_bases[k] = [(read_pos, score/255.0) for (read_pos, score) in this_read.modified_bases.get(('N', 0, v), [])]
            read_mod_bases[this_read.query_name] = this_read_mod_bases
    return read_mod_bases

thresh_dist = 2000
bin_width = 10
num_bins = 2*thresh_dist//bin_width
R_BINS = np.linspace(-thresh_dist, thresh_dist, num_bins+1)
R_CENTERS = 0.5 * (R_BINS[:-1] + R_BINS[1:])

def get_single_read_corr_func(in_read_mod_bases, mod_primary, mod_secondary, r_bins, min_mod_primary=0, max_mod_primary=1.01):
    r_probs = [[] for i in range(len(r_bins)-1)]

    primary_pos = [pos_prob[0] for pos_prob in in_read_mod_bases[mod_primary]
                   if ((pos_prob[1]>=min_mod_primary) and (pos_prob[1]<max_mod_primary))]
    if len(primary_pos)==0:
        return r_probs
    secondary_pos = []
    secondary_prob = []
    for pos_prob in in_read_mod_bases[mod_secondary]:
        secondary_pos.append(pos_prob[0])
        secondary_prob.append(pos_prob[1])
    if len(secondary_pos)==0:
        return r_probs

    r_mat = np.array(primary_pos)[:, np.newaxis] - np.array(secondary_pos)[np.newaxis, :]
    p_mat = np.tile(secondary_prob, (r_mat.shape[0], 1))

    for i_start in range(len(r_bins)-1):
        i_stop = i_start + 1
        r_start = r_bins[i_start]
        r_stop = r_bins[i_stop]
        mask_mat = (r_mat>=r_start) * (r_mat<r_stop)
        r_probs[i_start] = p_mat[mask_mat]
    return r_probs


def get_acc_corr_func(in_read_mod_bases, primary, secondary, min_mod, max_mod):
    acc_read_corr_func = []
    for this_read_mod_bases in tqdm(in_read_mod_bases.values()):
        acc_read_corr_func.append(get_single_read_corr_func(
            this_read_mod_bases, primary, secondary, R_BINS, min_mod_primary=min_mod, max_mod_primary=max_mod
        ))

    r_all_p = [[] for i_start in range(len(R_BINS)-1)]
    for l in acc_read_corr_func:
        for i in range(len(R_BINS)-1):
            r_all_p[i].extend(l[i])
    r_avg_p = [np.mean(l) for l in r_all_p]
    return np.array(r_avg_p)

########################################################################################################################

# chrom = '11'
for chrom in range(12, 23):
    print(f'chrom {chrom}')
    workspace = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT/chr{chrom}'
    bam_file = os.path.join(workspace, 'mAFiA.reads.bam')

    img_out = f'/home/adrian/img_out/psi-mAFiA_cross_talk/chr{chrom}'
    os.makedirs(img_out, exist_ok=True)

    read_mod_bases = get_read_mod_bases(bam_file)

    ### primary psi ###
    r_avg_p_high_psi = get_acc_corr_func(read_mod_bases, 'psi', 'm6A', min_mod=0.5, max_mod=1.01)
    r_avg_p_low_psi = get_acc_corr_func(read_mod_bases, 'psi', 'm6A', min_mod=0.0, max_mod=0.5)
    r_avg_p_high_m6A = get_acc_corr_func(read_mod_bases, 'm6A', 'psi', min_mod=0.5, max_mod=1.01)
    r_avg_p_low_m6A = get_acc_corr_func(read_mod_bases, 'm6A', 'psi', min_mod=0.0, max_mod=0.5)

    plt.figure(figsize=(10, 10))
    plt.subplot(2, 1, 1)
    plt.plot(R_CENTERS, r_avg_p_high_psi, label='$P(\psi)\geq0.5$')
    plt.plot(R_CENTERS, r_avg_p_low_psi, label='$P(\psi)<0.5$')
    plt.axvline(x=0, c='r', linestyle='--', alpha=0.5)
    plt.legend(loc='upper right')
    plt.xlabel('r($m^6A$) - r($\psi$)')
    plt.ylabel(r'$\langle P(m^6A) \rangle$')
    plt.subplot(2, 1, 2)
    plt.plot(R_CENTERS, r_avg_p_high_m6A, label='$P(m^6A)\geq0.5$')
    plt.plot(R_CENTERS, r_avg_p_low_m6A, label='$P(m^6A)<0.5$')
    plt.axvline(x=0, c='r', linestyle='--', alpha=0.5)
    plt.legend(loc='upper right')
    plt.xlabel('r($\psi$) - r($m^6A$)')
    plt.ylabel(r'$\langle P(\psi) \rangle$')
    plt.suptitle(f'chr{chrom}')
    plt.savefig(os.path.join(img_out, 'r_avg_p.png'), bbox_inches='tight')
    plt.close('all')

    df_out = pd.DataFrame({
        'r_center': np.int64(R_CENTERS),
        'avg_p_m6A_high_psi': r_avg_p_high_psi,
        'avg_p_m6A_low_psi': r_avg_p_low_psi,
        'avg_p_psi_high_m6A': r_avg_p_high_m6A,
        'avg_p_psi_low_m6A': r_avg_p_low_m6A,
    })
    df_out.to_csv(os.path.join(img_out, 'r_avg_p.tsv'), sep='\t', index=False, float_format='%.4f')

### sum over chrs ######################################################################################################
all_r_avg_p_high_psi = []
all_r_avg_p_low_psi = []
all_r_avg_p_high_m6A = []
all_r_avg_p_low_m6A = []

chr_start = 1
chr_end = 22

for chr in list(range(chr_start, chr_end + 1)) + ['X']:
    df = pd.read_csv(os.path.join('/home/adrian/img_out/psi-mAFiA_cross_talk', f'chr{chr}', 'r_avg_p.tsv'), sep='\t')
    all_r_avg_p_high_psi.append(df['avg_p_m6A_high_psi'].values)
    all_r_avg_p_low_psi.append(df['avg_p_m6A_low_psi'].values)
    all_r_avg_p_high_m6A.append(df['avg_p_psi_high_m6A'].values)
    all_r_avg_p_low_m6A.append(df['avg_p_psi_low_m6A'].values)
all_r_avg_p_high_psi = np.vstack(all_r_avg_p_high_psi).mean(axis=0)
all_r_avg_p_low_psi = np.vstack(all_r_avg_p_low_psi).mean(axis=0)
all_r_avg_p_high_m6A = np.vstack(all_r_avg_p_high_m6A).mean(axis=0)
all_r_avg_p_low_m6A = np.vstack(all_r_avg_p_low_m6A).mean(axis=0)

plt.figure(figsize=(10, 10))
plt.subplot(2, 1, 1)
plt.plot(R_CENTERS, all_r_avg_p_high_psi, label='$P(\psi)\geq0.5$')
plt.plot(R_CENTERS, all_r_avg_p_low_psi, label='$P(\psi)<0.5$')
# plt.axvline(x=0, c='r', linestyle='--', alpha=0.5)
plt.legend(loc='upper right')
plt.xlabel('r($m^6A$) - r($\psi$)')
plt.ylabel(r'$\langle P(m^6A) \rangle$')
plt.subplot(2, 1, 2)
plt.plot(R_CENTERS, all_r_avg_p_high_m6A, label='$P(m^6A)\geq0.5$')
plt.plot(R_CENTERS, all_r_avg_p_low_m6A, label='$P(m^6A)<0.5$')
# plt.axvline(x=0, c='r', linestyle='--', alpha=0.5)
plt.legend(loc='upper right')
plt.xlabel('r($\psi$) - r($m^6A$)')
plt.ylabel(r'$\langle P(\psi) \rangle$')
plt.suptitle(f'HEK293 WT\nchr {chr_start}-{chr_end}, X averaged')

plt.savefig(os.path.join('/home/adrian/img_out/psi-mAFiA_cross_talk', 'global_corr_func.png'), bbox_inches='tight')
plt.close('all')

