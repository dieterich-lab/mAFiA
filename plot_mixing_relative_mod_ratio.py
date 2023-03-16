import os
HOME = os.path.expanduser('~')
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# parser = argparse.ArgumentParser()
# parser.add_argument('--df_file')
# args = parser.parse_args()
# df_file = args.df_file
# df_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_multipleMotifs_noNorm_g3.tsv'

base_dataset = 'HEK293T-WT-100-rep1'
base_df_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_{}.tsv'.format(base_dataset)

comp_datasets = [
    'HEK293T-WT-0-rep2',
    'HEK293T-WT-50-rep3'
]
comp_df_files = [
    '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_{}.tsv'.format(comp_dataset) for comp_dataset in comp_datasets
]

img_out = os.path.join(HOME, 'img_out/MAFIA/WT_KO_mixing')
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

df_base = pd.read_csv(base_df_file, sep='\t', index_col=0)
dfs_comp = [pd.read_csv(comp_df_file, sep='\t', index_col=0) for comp_df_file in comp_df_files]

P_VAL_THRESH = 1.0E-10
COV_THRESH = 50
motifs = ['GGACA', 'GGACC', 'AGACT']

plt.figure(figsize=(15, 5))
for subplot_ind, this_motif in enumerate(motifs):
    df_base_thresh = df_base[
            (df_base['num_test_features']>=COV_THRESH) *
            (df_base['Pvalue']<P_VAL_THRESH) *
            (df_base['motif']==this_motif)
        ]
    dfs_comp_thresh = [
        df_comp[
            (df_comp['num_test_features']>=COV_THRESH) *
            (df_comp['Pvalue']<P_VAL_THRESH) *
            (df_comp['motif']==this_motif)
            ] for df_comp in dfs_comp
    ]
    common_idx = list(set.intersection(set(df_base_thresh.index), *[set(this_df.index) for this_df in dfs_comp_thresh]))

    x_vals = df_base_thresh.loc[common_idx]['mod_ratio'].values

    plt.subplot(1, 3, subplot_ind+1)
    for (df, ds) in zip(dfs_comp_thresh, comp_datasets):
        y_vals = df.loc[common_idx]['mod_ratio'].values
        plt.plot(x_vals, y_vals, '.', label=' '.join(ds.split('-')[1:3]))
        # plt.plot(x_vals, y_vals, 'o', mfc='none', label=' '.join(ds.split('-')[1:3]))
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('WT mod. ratio', fontsize=15)
    if subplot_ind==0:
        plt.ylabel('Mixture mod. ratio', fontsize=15)
    plt.legend(loc='upper left', fontsize=10)
    plt.title(this_motif, fontsize=15)

# plt.title(r'{} sites with ONT coverage$\geq${}'.format(df_motif.shape[0], coverage_thresh) + '\nCorrelation {:.2f}'.format(corr), fontsize=15)
# plt.subplots_adjust(top=0.8)
# plt.suptitle(this_motif, fontsize=20)
# plt.savefig(os.path.join(img_out, 'corr_glori_modRatio_{}_pValThresh{:.2E}.png'.format(this_motif, P_VAL_THRESH)), bbox_inches='tight')
# plt.close()