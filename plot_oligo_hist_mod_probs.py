import os
HOME = os.path.expanduser('~')
from glob import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def get_norm_counts(in_df, sel_motif):
    df_motif = in_df[
        (in_df['ref_motif']==sel_motif)
        # & (df['pred_motif']==sel_motif)
    ]

    ### histogram of mod prob per read ###
    num_bins = 100
    bin_edges = np.linspace(0, 1, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_counts, _ = np.histogram(df_motif['mod_prob'], bins=bin_edges)
    norm_counts = bin_counts / np.sum(bin_counts)

    return norm_counts, bin_centers, df_motif

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
training_dataset = 'WUE_combined'
testing_datasets = [
    'ISA_run1_A',
    'ISA_run1_m6A'
]
ds_colors = {
    'ISA_run1_A' : 'b',
    'ISA_run1_m6A' : 'r'
}

img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('oligo_validation_{}'.format(training_dataset)))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

COV_THRESH = 10

dfs = {}
for ds in testing_datasets:
    # paths = glob(os.path.join(results_dir, 'res_{}_{}_modProbPerRead.tsv.merged'.format(ds, training_dataset)))
    path = os.path.join(results_dir, 'debug_{}.tsv'.format(ds, training_dataset))
    df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})
    dfs[ds] = df

# motifs = list(set.intersection(*[set(df['ref_motif'].unique()) for df in dfs.values()]))
motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']

### plots ###
num_rows = 2
num_cols = 3
fig_width = 16
fig_height = 10
fig_hist = plt.figure(figsize=(fig_width, fig_height))
axes_hist = fig_hist.subplots(num_rows, num_cols).flatten()
for ds, df in dfs.items():
    for subplot_ind, this_motif in enumerate(motifs):
        if this_motif not in df['ref_motif'].unique():
            continue
        ds_norm_counts, ds_bin_centers, ds_motif = get_norm_counts(df, this_motif)
        axes_hist[subplot_ind].step(ds_bin_centers, ds_norm_counts, color=ds_colors[ds], label='{}, {} bases'.format(ds.split('_')[-1], len(ds_motif)))

        axes_hist[subplot_ind].set_xlabel('Mod. Prob.', fontsize=15)
        axes_hist[subplot_ind].set_ylabel('Norm. frequency', fontsize=15)
        axes_hist[subplot_ind].set_title('{}'.format(this_motif), fontsize=20)
        axes_hist[subplot_ind].set_xlim([-0.05, 1.05])
        axes_hist[subplot_ind].set_ylim([0, 0.5])
        axes_hist[subplot_ind].legend(loc='upper left', fontsize=12)

fig_hist.tight_layout()
fig_hist.savefig(os.path.join(img_out, 'hist_oligo_modProbs.png'), bbox_inches='tight')
