import os
HOME = os.path.expanduser('~')
import pandas as pd
import numpy as np
from tqdm import tqdm
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
# FMT = 'pdf'
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

# train_dataset = 'ISA-WUE'
# results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'

# train_dataset = 'CHEUI'
# results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/CHEUI/'

train_dataset = 'm6Anet'
results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6Anet'

test_datasets = [
    'WUE_batch3_AB_BA',
    'WUE_batch3_ABBA',
]
ds_names = {
    'WUE_batch3_AB_BA' : 'AB | BA',
    'WUE_batch3_ABBA' : 'AB - BA',
}
target_pattern = 'M0S1'

thresh_q = 70
block_size = 13
block_center = 6

img_out = os.path.join(HOME, 'img_out/MAFIA', 'res_train_{}_test_WUE_batch3_AB+BA_minimap'.format(train_dataset))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### load mod. probs. or collect from scratch ###########################################################################
########################################################################################################################
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

def import_mAFiA_res(res_dir):
    tol_pos = 2
    ds_pos_prob = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        # path = os.path.join(res_dir, 'res_train_{}_test_{}_q{}.tsv'.format(train_dataset, test_ds, thresh_q))
        path = os.path.join(res_dir, 'res_train_{}_test_{}_minimap.tsv'.format(train_dataset, test_ds))
        df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})
        unique_reads = df['read_id'].unique()
        ds_pos_prob[test_ds] = {}
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
            sel_refPos = refPos[sel_indices]
            sel_modProb = modProb[sel_indices]

            ds_pos_prob[test_ds][this_read] = list(zip(sel_refPos, sel_modProb))

    return ds_pos_prob

def import_m6Anet_res(res_dir):
    ds_pos_prob = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        path = os.path.join(res_dir, test_ds, 'data.indiv_proba.csv')
        df = pd.read_csv(path, sep=',', dtype = {
            'read_index': np.int64,
            'transcript_position': np.int64,
            'probability_modified': np.float64
        })
        unique_reads = df['read_index'].unique()
        ds_pos_prob[test_ds] = {}
        for this_read in tqdm(unique_reads):
            sub_df = df[df['read_index']==this_read]
            if len(sub_df)<2:
                continue

            refPos = np.int64(sub_df['transcript_position'].values)
            modProb = sub_df['probability_modified'].values

            ds_pos_prob[test_ds][this_read] = list(zip(refPos, modProb))

    return ds_pos_prob

def import_CHEUI_res(res_dir):
    ds_pos_prob = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        path = os.path.join(res_dir, test_ds, 'read_level_m6A_predictions.txt')
        df = pd.read_csv(path, sep='\t', names=['site', 'mod_prob', 'replicate'], dtype={'mod_prob': np.float64})

        contig = []
        ref_pos = []
        ref_motif = []
        read_id = []
        for this_site in df['site'].values:
            this_contig, this_start_pos, this_9mer, this_read_id = this_site.split('_')
            contig.append(this_contig)
            ref_pos.append(int(this_start_pos)+4)
            ref_motif.append(this_9mer[2:7])
            read_id.append(this_read_id)

        df['contig'] = contig
        df['ref_pos'] = ref_pos
        df['ref_motif'] = ref_motif
        df['read_id'] = read_id
        unique_reads = df['read_id'].unique()

        ds_pos_prob[test_ds] = {}
        for this_read in tqdm(unique_reads):
            sub_df = df[(df['read_id']==this_read) * (df['ref_motif']=='GGACT')]
            if len(sub_df) < 2:
                continue

            unique_contigs = sub_df['contig'].unique()
            longest_contig = unique_contigs[np.argmax([len(c) for c in unique_contigs])]
            sub_df = sub_df[sub_df['contig'] == longest_contig]

            refPos = np.int64(sub_df['ref_pos'].values)
            modProb = sub_df['mod_prob'].values
            ds_pos_prob[test_ds][this_read] = list(zip(refPos, modProb))

    return ds_pos_prob

if train_dataset=='ISA-WUE':
    ds_refPos_modProbs = import_mAFiA_res(results_dir)
elif train_dataset=='m6Anet':
    ds_refPos_modProbs = import_m6Anet_res(results_dir)
elif train_dataset=='CHEUI':
    ds_refPos_modProbs = import_CHEUI_res(results_dir)

########################################################################################################################
### construct probality mat ############################################################################################
########################################################################################################################
min_pos = 6
extent = 80
min_cycles = 2

for vmin in np.arange(0.1, 1, 0.1):

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
    num_samples = 30

    # num_samples = 10
    # ytick_pos = [0, 4, 9]

    fig_width = 5*cm
    fig_height = 5*cm

    num_rows = 2
    num_cols = 1

    # cax_yticks = np.linspace(vmin, 1, 3)
    xticks = np.arange(min_pos, extent, block_size)
    ytick_pos = np.concatenate([[0], np.arange(num_samples//2-1, num_samples, num_samples//2)])
    yticks = [f'read {i+1}' for i in ytick_pos]

    fig_single_read = plt.figure(figsize=(fig_width, fig_height))
    for subplot_ind, test_ds in enumerate(test_datasets):
        full_mat = ds_mat_modProb[test_ds]
        # sel_rows = np.arange(num_samples)
        sel_rows = np.argpartition(np.max(np.abs(np.diff(full_mat[:, np.arange(min_pos, extent, block_size)], axis=1)), axis=1), -num_samples)[-num_samples:]
        display_mat = full_mat[sel_rows]
        ax = fig_single_read.add_subplot(num_rows, num_cols, subplot_ind+1)
        im = ax.imshow(display_mat, vmin=vmin, vmax=1, cmap='plasma')
        # ax.set_xticks(xticks)
        if subplot_ind==(num_rows-1):
            ax.set_xticks(xticks)
            ax.set_xlabel('Aligned Position')
        else:
            ax.set_xticks([])
        ax.set_yticks([])
        ax.set_yticks(ytick_pos, yticks)
        ax.set_xlim([-0.5, extent-0.5])
        # ax.set_title(ds_names[test_ds])

    fig_single_read.tight_layout(rect=[0.1, 0.1, 0.8, 0.9])
    fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
    cax = fig_single_read.add_axes([0.85, 0.3, 0.03, 0.4])
    plt.colorbar(im, cax=cax)
    # plt.yticks([0.8, 0.9, 1.0])
    plt.yticks(np.linspace(vmin, 1.0, 3))
    plt.ylabel('P(m6A)', rotation=-90, labelpad=10)

    fig_single_read.savefig(os.path.join(img_out, f'single_read_vmin{vmin:.1f}.{FMT}'), **fig_kwargs)
    plt.close('all')

    ########################################################################################################################
    ### bar plots ##########################################################################################################
    ########################################################################################################################
    thresh_modProb = vmin
    fig_barplot = plt.figure(figsize=(fig_width, fig_height))
    ds_contrast = {}
    for subplot_ind, test_ds in enumerate(test_datasets):
        x_pos = np.arange(min_pos, extent, block_size)
        prob_mat = ds_mat_modProb[test_ds][:, x_pos]
        avg_mod_ratio = np.sum(prob_mat>=thresh_modProb, axis=0) / np.sum(prob_mat>0, axis=0)
        ds_contrast[test_ds] = np.mean(avg_mod_ratio[[2, 4]] - avg_mod_ratio[[1, 3]])

        ax = fig_barplot.add_subplot(num_rows, num_cols, subplot_ind+1)
        ax.bar(x_pos, avg_mod_ratio, width=5, label=f'{ds_names[test_ds]}')
        ax.axhline(y=0.5, c='r', linestyle='--', alpha=0.5)
        if subplot_ind==(num_rows-1):
            ax.set_xticks(xticks)
            ax.set_xlabel('Aligned Position')
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
    full_contrast = ds_contrast['WUE_batch3_AB_BA'] - ds_contrast['WUE_batch3_ABBA']
    fig_barplot.savefig(os.path.join(img_out, f'barplot_vmin{vmin:.1f}_contrast{full_contrast:.2f}.{FMT}'), **fig_kwargs)
    plt.close('all')

    ########################################################################################################################
    ### block distances ####################################################################################################
    ########################################################################################################################
    ds_diff = {}
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

        ds_diff[test_ds] = np.sum((freq[np.isin(nt_dists, [26, 52])] - freq[np.isin(nt_dists, [13, 39])]))

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

    inter_ds_diff = ds_diff['WUE_batch3_AB_BA'] - ds_diff['WUE_batch3_ABBA']
    # high = ds_diff['WUE_batch3_AB_BA']
    # low = ds_diff['WUE_batch3_ABBA']

    fig_dist_freq.savefig(os.path.join(img_out, f'distFreq_vmin{vmin:.1f}_diff{inter_ds_diff:.2f}.{FMT}'), **fig_kwargs)
    plt.close('all')