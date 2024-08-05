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
# lower_limit = 50
# upper_limit = 150
polyA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA'
# below_bam_file = os.path.join(polyA_dir, f'SHAM_{day}_reads_polyA_len_below_{lower_limit}.bam')
# above_bam_file = os.path.join(polyA_dir, f'SHAM_{day}_reads_polyA_len_above_{upper_limit}.bam')
polyA_conditions = [
    '$< 50$',
    # '[50, 150)',
    '$\geq 150$'
]
condition_colors = [
    'skyblue',
    # 'royalblue',
    'midnightblue'
]
bam_files = [
    os.path.join(polyA_dir, f'SHAM_{day}_reads_polyA_len_below_50.bam'),
    # os.path.join(polyA_dir, f'SHAM_{day}_reads_polyA_len_50_150.bam'),
    os.path.join(polyA_dir, f'SHAM_{day}_reads_polyA_len_above_150.bam'),
]

img_out = '/home/adrian/img_out/polyA_len'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def calc_mod_level_per_read(in_bam_file, thresh_mod_prob=192.0, num_bases=5):
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


mod_level_per_read = [
    calc_mod_level_per_read(this_bam_file)
    for this_bam_file in bam_files
]

num_bins = 20
max_bin = 1.0

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)

    for this_mod_level_per_read, this_polyA_cond, this_color in zip(mod_level_per_read, polyA_conditions, condition_colors):
        this_mod_vals = [val.get(this_mod, np.nan) for val in this_mod_level_per_read.values()]

        this_hist, bin_edges = np.histogram(this_mod_vals, range=[0, max_bin], bins=num_bins)
        this_hist = this_hist / np.sum(this_hist)
        bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
        plt.plot(bin_centers, this_hist, c=this_color, label=this_polyA_cond)

    # plt.hist(hist_below, bins, alpha=0.5, label=f'p(A)<{lower_limit} bps')
    # plt.hist(hist_above, bins, alpha=0.5, label=f'p(A)>={upper_limit} bps')
    # plt.yscale('log')
    plt.xlim([-0.01, 1.01])
    plt.ylim([0, 0.1])
    plt.legend(title='poly(A) len, bps')
    plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $ per read', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
plt.suptitle(f'SHAM {day}', fontsize=15)
plt.savefig(os.path.join(img_out, f'hist_avg_mods_per_read_grouped_by_polyA_len_{day}.png'), bbox_inches='tight')
