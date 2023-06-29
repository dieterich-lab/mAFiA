import os
HOME = os.path.expanduser('~')
import re
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import random
random.seed(10)
from random import sample
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
# FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
train_dataset = 'ISA-WUE'
test_datasets = [
    'WUE_batch3_AB_BA',
    'WUE_batch3_ABBA',
]
ds_names = {
    'WUE_batch3_AB_BA' : 'AB | BA',
    'WUE_batch3_ABBA' : 'ABBA',
}
target_pattern = 'M0S1'

thresh_q = 70
block_size = 13
block_center = 6

# out_pickle = os.path.join(results_dir, 'ISA_run4_M4M5_modProbs_q{}.pkl'.format(thresh_q))

img_out = os.path.join(HOME, 'img_out/MAFIA', 'res_train_{}_test_WUE_batch3_ABBA'.format(train_dataset))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### load mod. probs. or collect from scratch ###########################################################################
########################################################################################################################
dfs = {}
for test_dataset in test_datasets:
    path = os.path.join(results_dir, 'res_train_{}_test_{}_q{}.tsv'.format(train_dataset, test_dataset, thresh_q))
    this_df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})
    dfs[test_dataset] = this_df
    # dfs[test_dataset] = this_df[this_df['ref_motif']==this_df['pred_motif']]

def get_longest_continuous_indices(indices):
    grouped_indices = []
    this_group = [indices.pop(0)]
    while len(indices)>0:
        x = indices.pop(0)
        if x > (this_group[-1]+1):
            grouped_indices.append(this_group)
            this_group = [x]
        else:
            this_group.append(x)
    grouped_indices.append(this_group)
    return grouped_indices[np.argmax([len(l) for l in grouped_indices])]

### collect mod probs ###
tol_pos = 2
ds_refPos_modProbs = {}
for test_ds in test_datasets:
    print('\nDataset {}'.format(test_ds))
    df = dfs[test_ds]
    unique_reads = df['read_id'].unique()
    ds_refPos_modProbs[test_ds] = {}
    for this_read in tqdm(unique_reads):
        sub_df = df[df['read_id']==this_read]
        if len(sub_df)<2:
            continue

        readPos = np.int64(sub_df['read_pos'].values)
        refPos = np.int64(sub_df['ref_pos'].values)
        modProb = sub_df['mod_prob'].values

        diff_readPos = np.diff(readPos)
        diff_refPos = np.diff(refPos)
        if np.any(np.abs(diff_readPos-block_size)>tol_pos) or np.any(np.abs(diff_refPos-block_size)>tol_pos):
            continue
        sel_indices = np.where(np.abs(diff_readPos - diff_refPos)<=tol_pos)[0]
        if (len(sel_indices)==0):
            continue

        if not np.all(np.diff(sel_indices)==1):
            sel_indices = get_longest_continuous_indices(sel_indices)

        sel_indices = np.concatenate([sel_indices, [sel_indices[-1]+1]])
        # filtered_indices = np.concatenate([sel_indices, (sel_indices + 1)])
        sel_refPos = refPos[sel_indices]
        sel_modProb = modProb[sel_indices]

        ds_refPos_modProbs[test_ds][this_read] = list(zip(sel_refPos, sel_modProb))

########################################################################################################################
### construct probality mat ############################################################################################
########################################################################################################################
min_pos = 6
extent = 80
vmin = 0.8
min_cycles = 2

ds_mat_modProb = {}
for test_ds in test_datasets:
    refPos_modProbs = ds_refPos_modProbs[test_ds]
    # mod_probs = [pair for pair in ds_pattern_modProbs[test_ds][c_pattern] if np.max(pair)>=vmin]
    # print(test_ds, c_pattern, len(mod_probs), np.mean(np.vstack(mod_probs), axis=0))

    sample_refPos_modProbs = list(refPos_modProbs.values())
    mat_mod_prob = np.zeros([len(refPos_modProbs), extent])
    i = 0
    for pos_probs in sample_refPos_modProbs:
        refPos = np.array([tup[0] for tup in pos_probs])
        modProbs = np.array([tup[1] for tup in pos_probs])
        if np.sum(modProbs>vmin)<min_cycles:
            continue

        first_high_ind = np.where(modProbs>vmin)[0][0]
        # if first_high_ind==0:
        #     filtered_refPos = refPos
        #     filtered_modProbs = modProbs
        #     shifted_refPos = filtered_refPos - np.min(filtered_refPos) + min_pos + block_size
        # else:
        filtered_refPos = refPos[first_high_ind:]
        filtered_modProbs = modProbs[first_high_ind:]
        shifted_refPos = filtered_refPos - np.min(filtered_refPos) + min_pos

        extent_mask = shifted_refPos<extent
        shifted_refPos = shifted_refPos[extent_mask]
        filtered_modProbs = filtered_modProbs[extent_mask]

        # if np.max(shifted_pos)>=extent:
        #     continue
        # mat_mod_prob[i, filtered_pos] = filtered_probs
        mat_mod_prob[i, shifted_refPos] = filtered_modProbs
        i+=1

    mat_mod_prob = mat_mod_prob[~np.all(mat_mod_prob==0, axis=1)]
    ds_mat_modProb[test_ds] = mat_mod_prob
    print(f"{test_ds}: {len(mat_mod_prob)} samples")

########################################################################################################################
### visualize ##########################################################################################################
########################################################################################################################
num_samples = 40

fig_width = 5*cm
fig_height = 5*cm

num_rows = 2
num_cols = 1

# cax_yticks = np.linspace(vmin, 1, 3)

ytick_pos = [0, 19, 39]
yticks = [f'read {i+1}' for i in ytick_pos]

xticks = np.arange(min_pos, extent, 2*block_size)

fig_single_read = plt.figure(figsize=(fig_width, fig_height))
for subplot_ind, test_ds in enumerate(test_datasets):
    display_mat = ds_mat_modProb[test_ds][:num_samples]
    # display_mat = ds_mat_modProb[test_ds][sample(range(ds_mat_modProb[test_ds].shape[0]), num_samples)]
    # display_mat = display_mat[np.sum(display_mat>vmin, axis=1)>=3][:num_samples]
    ax = fig_single_read.add_subplot(num_rows, num_cols, subplot_ind+1)
    im = ax.imshow(display_mat, vmin=vmin, vmax=1, cmap='plasma')
    # ax.set_xticks(xticks)
    if subplot_ind==(num_rows-1):
        ax.set_xticks(xticks)
        ax.set_xlabel('Aligned position')
    else:
        ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks(ytick_pos, yticks)
    ax.set_title(ds_names[test_ds])

fig_single_read.tight_layout(rect=[0.1, 0.1, 0.8, 0.9])
fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
cax = fig_single_read.add_axes([0.85, 0.3, 0.03, 0.4])
plt.colorbar(im, cax=cax)
plt.yticks([0.8, 0.9, 1.0])
plt.ylabel('P(m6A)', rotation=-90, labelpad=10)

fig_single_read.savefig(os.path.join(img_out, f'single_read_q{thresh_q}.{FMT}'), **fig_kwargs)

########################################################################################################################
### bar plots ##########################################################################################################
########################################################################################################################
thresh_modProb = 0.8
fig_barplot = plt.figure(figsize=(fig_width, fig_height))
for subplot_ind, test_ds in enumerate(test_datasets):
    x_pos = np.arange(min_pos, extent, block_size)
    prob_mat = ds_mat_modProb[test_ds][:, x_pos]
    avg_mod_ratio = np.sum(prob_mat>=thresh_modProb, axis=0) / np.sum(prob_mat>0, axis=0)

    ax = fig_barplot.add_subplot(num_rows, num_cols, subplot_ind+1)
    ax.bar(x_pos, avg_mod_ratio, width=5, label=ds_names[test_ds])
    ax.axhline(y=0.5, c='r', linestyle='--', alpha=0.5)
    if subplot_ind==(num_rows-1):
        ax.set_xticks(xticks)
        ax.set_xlabel('Aligned position')
    else:
        ax.set_xticks([])
    ax.set_ylim([0.0, 1.05])
    ax.set_ylabel(f'Mod. Ratio')
    # if subplot_ind==1:
    #     ax.set_ylabel(f'Mod. Ratio')
        # ax.set_ylabel(f'% NTs with P(m6A)>={thresh_modProb:.1f}')
        # ax.set_ylabel('Mean P(m6A)')
    ax.legend(loc='upper right')
    # ax.set_yticks([])
    # ax.set_yticks(ytick_pos, yticks)
    # ax.set_title(ds_names[test_ds])
fig_barplot.savefig(os.path.join(img_out, f'barplot_q{thresh_q}.{FMT}'), **fig_kwargs)

########################################################################################################################
### violin plots #######################################################################################################
########################################################################################################################
# x_pos = np.arange(min_pos, extent, block_size)
#
# fig_vlnplot = plt.figure(figsize=(fig_width, fig_height))
# for row_ind, test_ds in enumerate(test_datasets):
#     vlnplot_mat = ds_mat_modProb[test_ds]
#     ax = fig_vlnplot.add_subplot(num_rows, num_cols, row_ind+1)
#     im = ax.violinplot(vlnplot_mat[:, x_pos])
#     if row_ind==(num_rows-1):
#         ax.set_xticks(xticks)
#         # ax.set_xlabel('Aligned pattern', fontsize=15)
#     else:
#         ax.set_xticks([])
#     # ax.set_ylim([0.0, 1.0])
#     # ax.set_yticks([])
#     # ax.set_yticks(ytick_pos, yticks)
#     ax.set_title(ds_names[test_ds])
# fig_vlnplot.savefig(os.path.join(img_out, f'vlnplot_q{thresh_q}.{FMT}'), **fig_kwargs)

########################################################################################################################
### block distances ####################################################################################################
########################################################################################################################
fig_dist_freq = plt.figure(figsize=(fig_width, fig_height))
for subplot_ind, test_ds in enumerate(test_datasets):
    mat_thresh = (ds_mat_modProb[test_ds]>=thresh_modProb)
    i_ind, j_ind = np.where(mat_thresh)
    dists = []
    for i in np.unique(i_ind):
        js = j_ind[i_ind == i]
        dists.extend(list(np.diff(js)))
    dist_counts = Counter(dists).most_common()
    # block_dists = [int(tup[0]/block_size) for tup in dist_counts]
    nt_dists = [int(tup[0]) for tup in dist_counts]
    counts = [tup[1] for tup in dist_counts]
    sorted_counts = np.array(counts)[np.argsort(nt_dists)]
    freq = sorted_counts / np.sum(sorted_counts)
    # block_dists = np.array(block_dists)[np.argsort(block_dists)]
    nt_dists = np.array(nt_dists)[np.argsort(nt_dists)]

    ax = fig_dist_freq.add_subplot(num_rows, num_cols, subplot_ind+1)
    ax.bar(np.arange(len(nt_dists)), freq, width=0.5, label=ds_names[test_ds])
    if subplot_ind==(num_rows-1):
        ax.set_xticks(np.arange(len(nt_dists)), nt_dists)
        ax.set_xlabel('m6A Distance')
    else:
        ax.set_xticks([])
    ax.set_ylabel('Norm. Frequency')
    ax.set_xlim([-0.5, 3.5])
    ax.set_ylim([0.0, 0.65])
    ax.set_yticks(np.arange(0, 0.65, 0.2))
    ax.legend(loc='upper right')
    # ax.set_title(ds_names[test_ds])

fig_dist_freq.savefig(os.path.join(img_out, f'distFreq_q{thresh_q}.{FMT}'), **fig_kwargs)