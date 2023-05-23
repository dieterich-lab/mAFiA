import os
from os.path import expanduser
HOME = expanduser('~')
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

thresh_n_reads = 50
thresh_pval = 1E-99

ds_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Christoph/m6anet/workflow_tx/inference'
out_dir = os.path.join(HOME, 'img_out/m6Anet')
ds_names = [
    '0_WT_100_IVT',
    # '25_WT_75_IVT',
    '50_WT_50_IVT',
    # '75_WT_25_IVT',
    '100_WT_0_IVT'
]

plt.figure(figsize=(15, 4))
for subplot_ind, this_ds_name in enumerate(ds_names):
    this_df = pd.read_csv(os.path.join(ds_dir, this_ds_name, 'data.site_proba_glori_filtered.csv'))
    thresh_df = this_df[
        (this_df['n_reads']>=thresh_n_reads)
        * (this_df['Pvalue']<thresh_pval)
        ]

    plt.subplot(1, 5, subplot_ind+1)
    plt.plot(thresh_df['GLORI'], thresh_df['mod_ratio'], '.')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('GLORI', fontsize=15)
    plt.ylabel('m6Anet', fontsize=15)
    plt.title('{}'.format(this_ds_name), fontsize=20)
plt.subplots_adjust(top=0.9)
plt.suptitle('m6Anet on HEK293\nn_reads$\\geq${}, p_val$\\leq${}'.format(thresh_n_reads, thresh_pval), fontsize=20)
plt.savefig(os.path.join(out_dir, 'all_ds_nreads{}_pval{}.png'.format(thresh_n_reads, thresh_pval)))