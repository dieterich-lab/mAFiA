import os
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

cluster_thresh = 0.90
# mafia_file = '/home/adrian/inference/rRNA/NR_003286_RNA18SN5_outlier_ratios_None_sigma0.97.tsv'
mafia_file = "/home/adrian/inference/rRNA/NR_003286_RNA18SN5_outlier_ratios_[all]_sigma{:.2f}.tsv".format(cluster_thresh)

base_file = '/home/adrian/Data/rRNA/Taoka_table.xlsx'
img_out = '/home/adrian/img_out'

df_mafia = pd.read_csv(mafia_file, sep='\t')

df_base = pd.read_excel(base_file, skiprows=[0])
df_base = df_base[df_base['rRNA']=='18S']

mod_counts = Counter(df_mafia['mod']).most_common()
mod_counts = [x for x in mod_counts if x[0]!='ND']

num_rows = 3
num_cols = 4
plt.figure(figsize=(num_cols*3, num_rows*3))
for subplot_ind, (this_mod, _) in enumerate(mod_counts):
    df_mafia_mod = df_mafia[df_mafia['mod']==this_mod]
    mod_ratios = []
    for pos in df_mafia_mod['stop']:
        r_base = df_base[df_base['Position assigned in this study']==pos]['Modification ratio'].values[0]
        r_mafia = df_mafia_mod[df_mafia_mod['stop']==pos]['ratio_outlier'].values[0] * 100
        mod_ratios.append([r_base, r_mafia])
    mod_ratios = np.vstack(mod_ratios)

    plt.subplot(num_rows, num_cols, subplot_ind+1)
    plt.plot(mod_ratios[:, 0], mod_ratios[:, 1], 'o', mfc='none', label=this_mod)
    plt.xlim([-5, 105])
    plt.ylim([-5, 105])
    plt.xticks(np.arange(0, 101, 25), np.arange(0, 101, 25))
    plt.yticks(np.arange(0, 101, 25), np.arange(0, 101, 25))
    if subplot_ind>=(num_rows-1)*num_cols:
        plt.xlabel('taoka', fontsize=15)
    if subplot_ind%num_cols==0:
        plt.ylabel('clustering', fontsize=15)
    plt.legend(loc='upper left', fontsize=12)
plt.suptitle('rRNA 18S\nClustering Threshold {:.2f}'.format(cluster_thresh), fontsize=20)
plt.savefig(os.path.join(img_out, 'mod_ratios_clustering_vs_taoka_thresh{:.2f}.png'.format(cluster_thresh)), bbox_inches='tight')