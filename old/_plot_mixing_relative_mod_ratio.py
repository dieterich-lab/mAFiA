import os
HOME = os.path.expanduser('~')
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# datasets = [
#     'HEK293A_WT',
#     'HEK293_IVT',
# ]
datasets = [
    #'HEK293T-WT-0-rep2',
    'HEK293T-Mettl3-KO-rep3',
    'HEK293T-WT-50-rep2',
    # 'HEK293T-WT-50-rep3',
    #'HEK293T-WT-100-rep1'
    'HEK293T-WT-rep3'
]
ds_colors = [
    'b',
    'y',
    'r'
]

mod_thresh = 0.95
dataset_names = datasets
# dataset_names = [' '.join(ds.split('-')[1:3]) for ds in datasets]
df_files = [
    '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_{}_modProbThresh{:.2f}.tsv'.format(dataset, mod_thresh) for dataset in datasets
]

img_out = os.path.join(HOME, 'img_out/MAFIA/Mettl3-mixing')
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

dfs = [pd.read_csv(df_file, sep='\t', index_col=0) for df_file in df_files]

P_VAL_THRESH = 1.0E-99
COV_THRESH = 100
motifs = ['GGACA', 'GGACC', 'AGACT']

plt.figure(figsize=(15, 5))
for subplot_ind, this_motif in enumerate(motifs):
    dfs_thresh = [
        df[
            (df['num_test_features']>=COV_THRESH) *
            (df['Pvalue']<P_VAL_THRESH) *
            (df['motif']==this_motif)
            ] for df in dfs
    ]
    common_idx = list(set.intersection(*[set(this_df.index) for this_df in dfs_thresh]))
    # print('{}: {} sites in common'.format(this_motif, len(common_idx)))

    x_vals = dfs_thresh[0].loc[common_idx]['Ratio'].values

    plt.subplot(1, 3, subplot_ind+1)
    for (df, ds, ds_name, ds_color) in zip(dfs_thresh, datasets, dataset_names, ds_colors):
        y_vals = df.loc[common_idx]['mod_ratio'].values
        corr = np.corrcoef(x_vals, y_vals)[0, 1]
        plt.plot(x_vals, y_vals, '.', color=ds_color, label='{}, corr. {:.2f}'.format(ds_name, corr))
        plt.plot(np.arange(0, 1.01, 0.1), np.arange(0, 1.01, 0.1), 'k--')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('GLORI mod. ratio', fontsize=15)
    if subplot_ind==0:
        plt.ylabel('Mixture mod. ratio', fontsize=15)
    plt.legend(loc='upper left', fontsize=10)
    plt.title('{}\n{} common sites'.format(this_motif, len(common_idx)), fontsize=15)

plt.subplots_adjust(top=0.8)
plt.suptitle('HEK293T WT / KO, Coverage $\geq$ {}'.format(COV_THRESH), fontsize=20)
plt.savefig(os.path.join(img_out, 'glori_modRatio_modProbThresh{:.2f}_pValThresh{:.2E}_covTHRESH{}.png'.format(mod_thresh, P_VAL_THRESH, COV_THRESH)), bbox_inches='tight')
plt.close()