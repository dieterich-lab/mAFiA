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


def calc_mod_level_per_read(in_bam_file, num_bases=5):
    read_avg_mod_level = {}
    with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
        for this_read in tqdm(in_bam.fetch()):
            read_avg_mod_level[this_read.query_name] = {}
            for this_mod in dict_mod_code.keys():
                mod_bases = this_read.modified_bases.get(('N', 0, dict_mod_code[this_mod]), [])

                if len(mod_bases)>=num_bases:
                    read_avg_mod_level[this_read.query_name][this_mod] = np.mean(np.partition([tup[1] for tup in mod_bases], -num_bases)[-num_bases:]) / 255.0
                else:
                    read_avg_mod_level[this_read.query_name][this_mod] = np.mean([tup[1] for tup in mod_bases]) / 255.0

    return read_avg_mod_level


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

########################################################################################################################
### define datasets ####################################################################################################
########################################################################################################################
ds = 'TAC'
conditions = ['SHAM_day56', 'TAC_day56']

# ds = 'HFpEF'
# conditions = ['ctrl', 'HFpEF']

# ds = 'Diet'
# conditions = ['WT_CD', 'WT_WD']
# # # conditions = ['M3KO_CD', 'M3KO_WD']
# # # conditions = ['WT_CD', 'M3KO_CD']
# # conditions = ['WT_WD', 'M3KO_WD']

# ds = 'CM'
# conditions = ['WT', 'M3KO']

polyA_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/polyA'
polyA_files = {this_cond: os.path.join(polyA_dir, f'{this_cond}_read2polyA_length.txt') for this_cond in conditions}

bam_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}'
if ds == 'TAC':
    bam_files = {this_cond: os.path.join(bam_dir, f'{this_cond}/chrALL.mAFiA.reads.bam') for this_cond in conditions}
else:
    bam_files = {this_cond: os.path.join(bam_dir, f'{this_cond}_merged/chrALL.mAFiA.reads.bam') for this_cond in conditions}

img_out = '/home/adrian/img_out/polyA_len'
os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### overall polyA distribution #########################################################################################
########################################################################################################################
df_polyA = {this_cond: pd.read_csv(polyA_files[this_cond], sep='\t', names=['read_id', 'polyA_len'])
            for this_cond in conditions}
dict_polyA = {this_cond: {k: v for k, v in df_polyA[this_cond][['read_id', 'polyA_len']].values}
              for this_cond in conditions}

plt.figure(figsize=(8, 4.5))
for cond_ind, this_cond in enumerate(conditions):
    plt.hist(dict_polyA[this_cond].values(), range=[0, 500], bins=50, label=this_cond, alpha=0.5, density=True)
plt.legend(fontsize=10)
plt.xlabel('polyA length (bps)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.title(ds, fontsize=15)
plt.savefig(os.path.join(img_out, f'hist_{ds}_{conditions[0]}_vs_{conditions[1]}_polyA_length.png'), bbox_inches='tight')

########################################################################################################################
### subset BAM files ###################################################################################################
########################################################################################################################
upper_lim = 150
lower_lim = 50

for this_cond, this_df in df_polyA.items():
    above_read_ids = this_df[this_df['polyA_len'] >= upper_lim]['read_id'].tolist()
    above_read_ids_file = os.path.join(polyA_dir, f'{this_cond}_read_ids_polyA_len_above_{upper_lim}.txt')
    above_bam_file = os.path.join(polyA_dir, f'{this_cond}_polyA_len_above_{upper_lim}.bam')

    with open(above_read_ids_file, 'w') as f_out:
        for this_read_id in above_read_ids:
            f_out.writelines(this_read_id + '\n')
    pysam.view("-N", above_read_ids_file, "-o", above_bam_file, bam_files[this_cond], catch_stdout=False)
    pysam.index(above_bam_file)

    below_read_ids = this_df[this_df['polyA_len'] < lower_lim]['read_id'].tolist()
    below_read_ids_file = os.path.join(polyA_dir, f'{this_cond}_read_ids_polyA_len_below_{lower_lim}.txt')
    below_bam_file = os.path.join(polyA_dir, f'{this_cond}_polyA_len_below_{lower_lim}.bam')
    with open(below_read_ids_file, 'w') as f_out:
        for this_read_id in below_read_ids:
            f_out.writelines(this_read_id + '\n')
    pysam.view("-N", below_read_ids_file, "-o", below_bam_file, bam_files[this_cond], catch_stdout=False)
    pysam.index(below_bam_file)

########################################################################################################################
### top mod. levels per read ###########################################################################################
########################################################################################################################
mod_level_per_read = {
    this_cond: {
        polyA_suffix: calc_mod_level_per_read(os.path.join(polyA_dir, f'{this_cond}_polyA_len_{polyA_suffix}.bam'))
        for polyA_suffix in [f'below_{lower_lim}', f'above_{upper_lim}']
    }
    for this_cond in conditions
}

num_bins = 20
max_bin = 1.0

cond_ls = {
    this_cond: this_ls for this_cond, this_ls in zip(conditions, ['-', '--'])
}

polyA_text = [
    f'$< {lower_lim}$',
    # '[50, 150)',
    f'$\geq {upper_lim}$'
]
polyA_colors = [
    'skyblue',
    # 'royalblue',
    'midnightblue'
]

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    for this_cond in conditions:
        for this_mod_level_per_read, this_polyA_text, this_polyA_color in \
                zip(mod_level_per_read[this_cond].values(), polyA_text, polyA_colors):
            this_mod_vals = [val.get(this_mod, np.nan) for val in this_mod_level_per_read.values()]
            this_hist, bin_edges = np.histogram(this_mod_vals, range=[0, max_bin], bins=num_bins)
            this_hist = this_hist / np.sum(this_hist)
            bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
            plt.plot(bin_centers, this_hist, c=this_polyA_color, ls=cond_ls[this_cond], label=f'{this_cond}, {this_polyA_text}')

    plt.xlim([-0.01, 1.01])
    plt.ylim([0, 0.1])
    plt.legend(title='poly(A) len, bps')
    plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $ per read', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
plt.suptitle(f'{conditions[0]} vs {conditions[1]}', fontsize=15)
plt.savefig(os.path.join(img_out, f'mod_per_read_{ds}_{conditions[0]}_vs_{conditions[1]}.png'), bbox_inches='tight')
