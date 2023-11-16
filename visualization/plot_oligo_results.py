import os
HOME = os.path.expanduser('~')
import pandas as pd
import numpy as np
import matplotlib as mpl
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
random.seed(0)
from random import sample
import re

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

xtick_spacing = 0.5
ytick_spacing = 0.1


def get_df_samples(in_df, sample_size=10000):
    return in_df.iloc[sample(range(len(in_df)), sample_size)]

def get_norm_counts(in_df, sel_motif):
    df_motif_all = in_df[
        (in_df['ref_motif']==sel_motif)
        # & (df['pred_motif']==sel_motif)
    ]

    df_motif = get_df_samples(df_motif_all)

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

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/results'
# training_dataset = 'ISA'
# training_dataset = 'WUE'
training_dataset = 'ISA-WUE'

# testing_datasets = ['ISA_A', 'ISA_m6A']
# testing_datasets = ['WUE_A', 'WUE_m6A']
testing_datasets = ['ISA-WUE_A', 'ISA-WUE_m6A']

ds_colors = {ds : c for ds, c in zip(testing_datasets, ['b', 'r'])}

img_out = os.path.join(HOME, 'img_out/NCOMMS_rev', 'res_train_{}_test_{}'.format(training_dataset, training_dataset))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

dfs = {}
for ds in testing_datasets:
    path = os.path.join(results_dir, 'res_train_{}_test_{}_q70.tsv'.format(training_dataset, ds))
    # path = os.path.join(results_dir, 'debug_res_train_{}_test_{}.tsv'.format(training_dataset, ds))
    df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})
    dfs[ds] = df

# motifs = list(set.intersection(*[set(df['ref_motif'].unique()) for df in dfs.values()]))
motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']

### plots ###
num_rows = 2
num_cols = 3
fig_width = 8*cm
fig_height = fig_width / gr

x_max = 1.0
xticks = np.round(np.linspace(0, x_max, 3), 3)

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
        # this_motif_bin_A_m6A_counts['{}_counts'.format(A_m6A)] = ds_real_counts
        this_motif_bin_A_m6A_counts['{}_counts'.format(A_m6A)] = ds_norm_counts

        if A_m6A=='A':
            label = 'UNM'
        elif A_m6A=='m6A':
            label = 'MOD'

        # y_max = max(y_max, (ds_norm_counts.max() // 0.05 + 1) * 0.05)
        # y_max = 0.2
        # axes_hist[subplot_ind].step(ds_bin_centers, ds_norm_counts, color=ds_colors[ds], label='{}, {} NTs'.format(ds.split('_')[-1], len(ds_motif)))
        # axes_hist[subplot_ind].step(ds_bin_centers, ds_norm_counts, color=ds_colors[ds], label='{}'.format(ds.split('_')[-1]))
        axes_hist[subplot_ind].step(ds_bin_centers, ds_real_counts, color=ds_colors[ds], label=label)

    if subplot_ind>=(num_rows-1)*num_cols:
        # axes_hist[subplot_ind].set_xlabel('Mod. Prob.')
        axes_hist[subplot_ind].set_xticks(xticks)
    else:
        axes_hist[subplot_ind].set_xticks([])

    # if subplot_ind%num_cols==0:
    #     axes_hist[subplot_ind].set_ylabel('Norm. frequency')

    # yticks = np.round(np.linspace(0, y_max, 3), 3)
    # axes_hist[subplot_ind].set_yticks(yticks)

    axes_hist[subplot_ind].set_title('{}'.format(this_motif), pad=-10)
    axes_hist[subplot_ind].set_xlim([-0.01, 1.01])
    # axes_hist[subplot_ind].set_ylim([-0.001, y_max])
    if subplot_ind==0:
        axes_hist[subplot_ind].legend(loc='upper left')
    motif_bin_A_m6A_counts[this_motif] = this_motif_bin_A_m6A_counts
# fig_hist.tight_layout(rect=[0, 0.03, 1, 0.9])
# fig_hist.suptitle('Trained on {}, tested on {}'.format(training_dataset, '_'.join(testing_datasets[0].split('_')[:-1])), fontsize=25)
fig_hist.tight_layout(pad=0.5)
fig_hist.savefig(os.path.join(img_out, f'hist_oligo_modProbs.{FMT}'), **fig_kwargs)

########################################################################################################################
### precision-recall curve #############################################################################################
########################################################################################################################
all_aucs = {}
fig_prc = plt.figure(figsize=(fig_width, fig_height))
axes_prc = fig_prc.subplots(num_rows, num_cols).flatten()
for subplot_ind, this_motif in enumerate(motifs):
    this_motif_recall, this_motif_precision = get_precision_recall_curve(motif_bin_A_m6A_counts[this_motif])
    this_motif_auc = get_auc(this_motif_recall, this_motif_precision)
    all_aucs[this_motif] = this_motif_auc
    axes_prc[subplot_ind].plot(this_motif_recall, this_motif_precision, label='AUC {:.2f}'.format(this_motif_auc))
    # if subplot_ind >= (num_rows - 1) * num_cols:
    #     axes_prc[subplot_ind].set_xlabel('Recall')
    # if subplot_ind % num_cols == 0:
    #     axes_prc[subplot_ind].set_ylabel('precision')
    axes_prc[subplot_ind].set_title('{}'.format(this_motif), pad=-10)
    axes_prc[subplot_ind].set_xlim([-0.01, 1.01])
    axes_prc[subplot_ind].set_ylim([-0.01, 1.01])
    axes_prc[subplot_ind].legend(loc='lower left')

    x_max = 1.0
    y_max = 1.0
    xticks = np.round(np.linspace(0, x_max, 3), 3)
    yticks = np.round(np.linspace(0, y_max, 3), 3)
    if subplot_ind>=(num_rows-1)*num_cols:
        axes_prc[subplot_ind].set_xticks(xticks)
    else:
        axes_prc[subplot_ind].set_xticks([])
    # if subplot_ind % num_cols == 0:
    #     axes_prc[subplot_ind].set_yticks(yticks)
    # else:
    #     axes_prc[subplot_ind].set_yticks([])
    axes_prc[subplot_ind].set_yticks(yticks)

# fig_prc.tight_layout(rect=[0, 0.03, 1, 0.9])
# fig_prc.suptitle('Trained on {}, tested on {}\nMean AUC {:.2f}'.format(training_dataset, '_'.join(testing_datasets[0].split('_')[:-1]), np.mean(list(all_aucs.values()))), fontsize=25)
fig_prc.tight_layout(pad=0.5)
fig_prc.savefig(os.path.join(img_out, f'prc_oligo_modProbs.{FMT}'), **fig_kwargs)

# plt.close('all')

### output mean AUC ###
mean_auc = np.mean(list(all_aucs.values()))
with open(os.path.join(img_out, 'mean_auc.txt'), 'w') as out_file:
    out_file.write('{:.3f}'.format(mean_auc))

### visualize single-read ###
# num_samples = 50
# max_pos = 150
# first_pos = 10
# block_size = 21
# xticks = np.arange(first_pos, max_pos, block_size)
# vmin = 0.8
# cax_yticks = np.arange(vmin, 1.01, 0.1)
#
# fig_single_read = plt.figure(figsize=(16, 10))
# for subplot_ind, ds in enumerate(testing_datasets):
#     df = dfs[ds]
#     read_ids = df['read_id'].unique()
#     pos_modProb = {}
#     for read_id in read_ids:
#         df_read = df[df['read_id']==read_id]
#         # if (len(df_read)>=6) and (int(df_read['contig'].unique()[0].split('_')[0].lstrip('block'))==len(df_read)):
#         if (len(df_read) >= 6) and (len(re.findall('-', df_read['contig'].unique()[0])) == len(df_read)):
#             pos_modProb[read_id] = df_read[['ref_pos', 'mod_prob']].values
#
#     mat_mod_prob = np.zeros([num_samples, max_pos])
#     ids = []
#     for i, k in enumerate(sample(list(pos_modProb.keys()), num_samples)):
#         ids.append(k.split('-')[-1])
#         qp = pos_modProb[k]
#         qs = np.int32(qp[:, 0])
#         ps = qp[:, 1]
#         qs = qs - np.min(qs) + first_pos
#         ps = ps[qs<max_pos]
#         qs = qs[qs<max_pos]
#         mat_mod_prob[i, qs] = ps
#
#     ax = fig_single_read.add_subplot(2, 1, subplot_ind+1)
#     im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1)
#     ax.set_xticks(xticks, fontsize=8)
#     ax.set_yticks(np.arange(num_samples), ids, fontsize=8)
#     if subplot_ind==1:
#         ax.set_xlabel('Aligned pos (NTs)', fontsize=20)
#     ax.set_ylabel(' '.join(ds.split('_')), fontsize=25, rotation=-90, labelpad=30)
#     ax.yaxis.set_label_position('right')
# fig_single_read.tight_layout(rect=[0.1, 0.2, 0.8, 0.8])
# fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
# cax = fig_single_read.add_axes([0.85, 0.35, 0.02, 0.3])
# plt.colorbar(im, cax=cax)
# plt.ylabel('P(m6A)', fontsize=20, rotation=-90, labelpad=30)
# plt.yticks(cax_yticks, cax_yticks)
#
# fig_single_read.savefig(os.path.join(img_out, 'single_read_oligo_modProbs.eps'))

plt.close('all')