import os
HOME = os.path.expanduser('~')
from glob import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random
random.seed(0)
from random import sample


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

    return bin_counts, norm_counts, bin_centers, df_motif

def get_precision_recall_curve(bin_A_m6A_counts):
    tp = np.cumsum(bin_A_m6A_counts['m6A_counts'][::-1])[::-1]
    fa = np.cumsum(bin_A_m6A_counts['A_counts'][::-1])[::-1]
    recall = tp / tp.max()
    precision = tp / (tp + fa)
    recall = np.concatenate([recall, [0.0]])
    precision = np.concatenate([precision, [1.0]])
    return recall, precision

def get_auc(rec, prec):
    return np.sum((prec[1:] + prec[:-1]) * 0.5 * np.abs(np.diff(rec)))

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
training_dataset = 'WUE_batches1+2'
testing_datasets = [
    'ISA_run1_A',
    'ISA_run1_m6A'
]
ds_colors = {
    'ISA_run1_A' : 'b',
    'ISA_run1_m6A' : 'r'
}

img_out = os.path.join(HOME, 'img_out/MAFIA', 'oligo_train_{}_test_{}'.format(training_dataset, '_'.join(testing_datasets[0].split('_')[:-1])))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

dfs = {}
for ds in testing_datasets:
    path = os.path.join(results_dir, 'res_train_{}_test_{}.tsv'.format(training_dataset, ds))
    # path = os.path.join(results_dir, 'debug_res_train_{}_test_{}.tsv'.format(training_dataset, ds))
    df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})
    dfs[ds] = df

# motifs = list(set.intersection(*[set(df['ref_motif'].unique()) for df in dfs.values()]))
motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']

### plots ###
num_rows = 2
num_cols = 3
fig_width = 16
fig_height = 12

### histogram of prob. distribution ###
fig_hist = plt.figure(figsize=(fig_width, fig_height))
axes_hist = fig_hist.subplots(num_rows, num_cols).flatten()
motif_bin_A_m6A_counts = {}
for subplot_ind, this_motif in enumerate(motifs):
    this_motif_bin_A_m6A_counts = {}
    y_max = 0
    for ds, df in dfs.items():
        if this_motif not in df['ref_motif'].unique():
            continue
        ds_real_counts, ds_norm_counts, ds_bin_centers, ds_motif = get_norm_counts(df, this_motif)
        this_motif_bin_A_m6A_counts['bin_centers'] = ds_bin_centers
        A_m6A = ds.split('_')[-1]
        this_motif_bin_A_m6A_counts['{}_counts'.format(A_m6A)] = ds_real_counts

        y_max = max(y_max, (ds_norm_counts.max() // 0.05 + 1) * 0.05)
        axes_hist[subplot_ind].step(ds_bin_centers, ds_norm_counts, color=ds_colors[ds], label='{}, {} NTs'.format(ds.split('_')[-1], len(ds_motif)))
    if subplot_ind>=(num_rows-1)*num_cols:
        axes_hist[subplot_ind].set_xlabel('Mod. Prob.', fontsize=15)
    if subplot_ind%num_cols==0:
        axes_hist[subplot_ind].set_ylabel('Norm. frequency', fontsize=15)
    axes_hist[subplot_ind].set_title('{}'.format(this_motif), fontsize=20)
    axes_hist[subplot_ind].set_xlim([-0.05, 1.05])
    axes_hist[subplot_ind].set_ylim([0, y_max])
    axes_hist[subplot_ind].legend(loc='upper left', fontsize=12)
    motif_bin_A_m6A_counts[this_motif] = this_motif_bin_A_m6A_counts
fig_hist.tight_layout(rect=[0, 0.03, 1, 0.9])
fig_hist.suptitle('Trained on {}\nTested on {}'.format(training_dataset, '_'.join(testing_datasets[0].split('_')[:-1])), fontsize=25)
fig_hist.savefig(os.path.join(img_out, 'hist_oligo_modProbs.png'), bbox_inches='tight')

### precision-recall curve ###
fig_prc = plt.figure(figsize=(fig_width, fig_height))
axes_prc = fig_prc.subplots(num_rows, num_cols).flatten()
for subplot_ind, this_motif in enumerate(motifs):
    this_motif_recall, this_motif_precision = get_precision_recall_curve(motif_bin_A_m6A_counts[this_motif])
    this_motif_auc = get_auc(this_motif_recall, this_motif_precision)
    axes_prc[subplot_ind].plot(this_motif_recall, this_motif_precision, label='AUC {:.2f}'.format(this_motif_auc))
    if subplot_ind >= (num_rows - 1) * num_cols:
        axes_prc[subplot_ind].set_xlabel('Recall', fontsize=15)
    if subplot_ind % num_cols == 0:
        axes_prc[subplot_ind].set_ylabel('precision', fontsize=15)
    axes_prc[subplot_ind].set_title('{}'.format(this_motif), fontsize=20)
    axes_prc[subplot_ind].set_xlim([-0.05, 1.05])
    axes_prc[subplot_ind].set_ylim([-0.05, 1.05])
    axes_prc[subplot_ind].legend(loc='lower left', fontsize=20)
fig_prc.tight_layout(rect=[0, 0.03, 1, 0.9])
fig_prc.suptitle('Trained on {}\nTested on {}'.format(training_dataset, '_'.join(testing_datasets[0].split('_')[:-1])), fontsize=25)
fig_prc.savefig(os.path.join(img_out, 'prc_oligo_modProbs.png'), bbox_inches='tight')

# plt.close('all')

### visualize single-read ###
num_samples = 50
max_pos = 150
first_pos = 10
block_size = 21
xticks = np.arange(first_pos, max_pos, block_size)
vmin = 0.8
cax_yticks = np.arange(vmin, 1.01, 0.1)

fig_single_read = plt.figure(figsize=(16, 10))
for subplot_ind, ds in enumerate(['ISA_run1_A', 'ISA_run1_m6A']):
    df = dfs[ds]
    read_ids = df['read_id'].unique()
    pos_modProb = {}
    for read_id in read_ids:
        df_read = df[df['read_id']==read_id]
        if (len(df_read)>=6) and (int(df_read['contig'].unique()[0].split('_')[0].lstrip('block'))==len(df_read)):
            pos_modProb[read_id] = df_read[['t_pos', 'mod_prob']].values

    mat_mod_prob = np.zeros([num_samples, max_pos])
    ids = []
    for i, k in enumerate(sample(list(pos_modProb.keys()), num_samples)):
        ids.append(k.split('-')[-1])
        qp = pos_modProb[k]
        qs = np.int32(qp[:, 0])
        ps = qp[:, 1]
        qs = qs - np.min(qs) + first_pos
        ps = ps[qs<max_pos]
        qs = qs[qs<max_pos]
        mat_mod_prob[i, qs] = ps

    ax = fig_single_read.add_subplot(2, 1, subplot_ind+1)
    im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1)
    ax.set_xticks(xticks, fontsize=8)
    ax.set_yticks(np.arange(num_samples), ids, fontsize=8)
    if subplot_ind==1:
        ax.set_xlabel('Aligned pos (NTs)', fontsize=20)
    ax.set_ylabel(' '.join(ds.split('_')), fontsize=25, rotation=-90, labelpad=30)
    ax.yaxis.set_label_position('right')
fig_single_read.tight_layout(rect=[0.1, 0.2, 0.8, 0.8])
fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
cax = fig_single_read.add_axes([0.85, 0.35, 0.02, 0.3])
plt.colorbar(im, cax=cax)
plt.ylabel('P(m6A)', fontsize=20, rotation=-90, labelpad=30)
plt.yticks(cax_yticks, cax_yticks)

fig_single_read.savefig(os.path.join(img_out, 'single_read_oligo_modProbs.png'))