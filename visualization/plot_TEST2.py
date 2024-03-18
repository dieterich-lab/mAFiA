import os
HOME = os.path.expanduser('~')
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.metrics import auc

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

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/TEST2'

# mAFiA_results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/results'
# mAFiA_results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/WUE_ABBA/mAFiA'
mAFiA_results_dir = os.path.join(source_data_dir, 'mAFiA')

# CHEUI_results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/CHEUI/oligo'
# CHEUI_results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/WUE_ABBA/CHEUI'
CHEUI_results_dir = os.path.join(source_data_dir, 'CHEUI')

# m6Anet_results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/m6Anet/oligo'
# m6Anet_results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/WUE_ABBA/m6Anet'
m6Anet_results_dir = os.path.join(source_data_dir, 'm6Anet')

test_datasets = [
    'WUE_batch3_AB_BA',
    # 'WUE_batch3_ABBA',
]
ds_names = {
    'WUE_batch3_AB_BA' : 'AB | BA',
    # 'WUE_batch3_ABBA' : 'AB - BA',
}
target_pattern = 'M0S1'

thresh_q = 70
block_size = 13
block_center = 6

# img_out = os.path.join(HOME, 'img_out/MAFIA', 'res_train_{}_test_WUE_batch3_AB+BA_minimap'.format(train_dataset))
# img_out = os.path.join(HOME, 'img_out/NCOMMS_rev', 'WUE_ABBA')
img_out = '/home/adrian/NCOMMS_revision/images/TEST2'
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


def filter_ref_pos_mod_prob(refPos, modProb, readPos=[], tol_pos=2):
    if len(readPos)==0:
        readPos = refPos
    diff_readPos = np.diff(readPos)
    diff_refPos = np.diff(refPos)
    # if np.any(np.abs(diff_readPos - block_size) > tol_pos) or np.any(np.abs(diff_refPos - block_size) > tol_pos):
    #     return [], []
    sel_indices = np.where(np.abs(diff_readPos - diff_refPos) <= tol_pos)[0]
    if (len(sel_indices) == 0):
        return [], []

    if not np.all(np.diff(sel_indices) == 1):
        sel_indices = get_longest_continuous_indices(list(sel_indices))

    sel_indices = np.concatenate([sel_indices, [sel_indices[-1] + 1]])
    sel_refPos = refPos[sel_indices]
    sel_modProb = modProb[sel_indices]

    return sel_refPos, sel_modProb


def import_mAFiA_res(res_dir):
    ds_pos_prob = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        # path = os.path.join(res_dir, 'res_train_{}_test_{}_q{}.tsv'.format(train_dataset, test_ds, thresh_q))
        # path = os.path.join(res_dir, 'res_train_ISA-WUE_test_{}_minimap.tsv'.format(test_ds))
        path = os.path.join(res_dir, 'pos_modProb.tsv')

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

            sel_refPos, sel_modProb = filter_ref_pos_mod_prob(refPos, modProb, readPos)
            # sel_refPos, sel_modProb = filter_ref_pos_mod_prob(refPos, modProb)
            if len(sel_refPos):
                ds_pos_prob[test_ds][this_read] = list(zip(sel_refPos, sel_modProb))

    return ds_pos_prob

def import_m6Anet_res(res_dir):
    ds_pos_prob = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        path = os.path.join(res_dir, 'data.indiv_proba.csv')
        df = pd.read_csv(path, sep=',', dtype = {
            'read_index': np.int64,
            'transcript_position': np.int64,
            'probability_modified': np.float64
        })

        df_read_index = pd.read_csv(os.path.join(res_dir, 'renata.fasta.index.fai'), sep='\t', names=['read_id'], usecols=[0])

        unique_read_indices = df['read_index'].unique()
        ds_pos_prob[test_ds] = {}
        for this_read_index in tqdm(unique_read_indices):
            sub_df = df[df['read_index']==this_read_index]
            if len(sub_df)<2:
                continue

            refPos = np.int64(sub_df['transcript_position'].values)
            modProb = sub_df['probability_modified'].values

            sel_refPos, sel_modProb = filter_ref_pos_mod_prob(refPos, modProb)
            if len(sel_refPos):
                this_read_id = df_read_index.loc[this_read_index].values[0]
                ds_pos_prob[test_ds][this_read_id] = list(zip(refPos, modProb))

    return ds_pos_prob

def import_CHEUI_res(res_dir):
    ds_pos_prob = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        path = os.path.join(res_dir, 'read_level_m6A_predictions.txt')
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

            # unique_contigs = sub_df['contig'].unique()
            # longest_contig = unique_contigs[np.argmax([len(c) for c in unique_contigs])]
            # sub_df = sub_df[sub_df['contig'] == longest_contig]

            refPos = np.int64(sub_df['ref_pos'].values)
            modProb = sub_df['mod_prob'].values

            sel_refPos, sel_modProb = filter_ref_pos_mod_prob(refPos, modProb)
            # sel_refPos, sel_modProb = refPos, modProb

            if len(sel_refPos):
                ds_pos_prob[test_ds][this_read] = list(zip(sel_refPos, sel_modProb))

    return ds_pos_prob


def get_even_odd_probs(read_mod_prob, block_size=13):
    ref_pos, mod_prob = np.vstack(read_mod_prob).T
    ref_pos = np.int32(ref_pos)
    highest_p_pos = ref_pos[np.argmax(mod_prob)]
    ref_pos - highest_p_pos
    even_mask = np.array([(pos / block_size) % 2 == 0 for pos in (ref_pos - highest_p_pos)])
    odd_mask = np.array([(pos / block_size) % 2 == 1 for pos in (ref_pos - highest_p_pos)])
    even_probs = mod_prob[even_mask]
    odd_probs = mod_prob[odd_mask]

    if len(even_probs)!=len(odd_probs):
        min_len = min(len(even_probs), len(odd_probs))
        return even_probs[:min_len], odd_probs[:min_len]

    return even_probs, odd_probs


def get_ds_even_odd_probs(ds_refPos_modProbs):
    even_odd_probs = {}
    even_odd_probs['even'] = []
    even_odd_probs['odd'] = []
    for read_refPos_modProbs in ds_refPos_modProbs.values():
        even_ps, odd_ps = get_even_odd_probs(read_refPos_modProbs)
        even_odd_probs['even'].extend(even_ps)
        even_odd_probs['odd'].extend(odd_ps)
    return even_odd_probs


def filter_by_readPos(in_modProbs, filter_modProbs):
    out_modProbs = {}
    for ds in in_modProbs.keys():
        out_modProbs[ds] = {}
        for readID in in_modProbs[ds].keys():
            if readID in filter_modProbs[ds].keys():
                read_list = []
                filter_pos = [x[0] for x in filter_modProbs[ds][readID]]
                for tup in in_modProbs[ds][readID]:
                    if tup[0] in filter_pos:
                        read_list.append(tup)
                if len(read_list):
                    out_modProbs[ds][readID] = read_list
        samples_before = np.sum([len(in_modProbs[ds][read]) for read in in_modProbs[ds].keys()])
        samples_after = np.sum([len(out_modProbs[ds][read]) for read in out_modProbs[ds].keys()])
        print(f'Before: {samples_before} samples')
        print(f'After: {samples_after}')
    return out_modProbs


mAFiA_refPos_modProbs = import_mAFiA_res(mAFiA_results_dir)
mAFiA_even_odd_probs = get_ds_even_odd_probs(mAFiA_refPos_modProbs['WUE_batch3_AB_BA'])
num_samples = len(mAFiA_even_odd_probs['even']) + len(mAFiA_even_odd_probs['odd'])

m6Anet_refPos_modProbs = import_m6Anet_res(m6Anet_results_dir)
# m6Anet_filtered_refPos_modProbs = filter_by_readPos(m6Anet_refPos_modProbs, mAFiA_refPos_modProbs)
m6Anet_filtered_refPos_modProbs = m6Anet_refPos_modProbs
m6Anet_even_odd_probs = get_ds_even_odd_probs(m6Anet_filtered_refPos_modProbs['WUE_batch3_AB_BA'])

CHEUI_refPos_modProbs = import_CHEUI_res(CHEUI_results_dir)
# CHEUI_filtered_refPos_modProbs = filter_by_readPos(CHEUI_refPos_modProbs, mAFiA_refPos_modProbs)
CHEUI_filtered_refPos_modProbs = CHEUI_refPos_modProbs
CHEUI_even_odd_probs = get_ds_even_odd_probs(CHEUI_filtered_refPos_modProbs['WUE_batch3_AB_BA'])

with open(os.path.join(img_out, 'num_samples_ABBA.tsv'), 'w') as h:
    h.write(f'GGACU\t{num_samples}\n')


########################################################################################################################
### site-level prob ####################################################################################################
########################################################################################################################
mAFiA_site_level_file = '/home/adrian/NCOMMS_revision/source_data/TEST2/mAFiA/mAFiA.sites.bed'
df_mAFiA_site_level = pd.read_csv(mAFiA_site_level_file, sep='\t')

x_start = 82
x_end = 151
x_ticks = np.arange(84, x_end, 13)

fig_site_level, ax = plt.subplots(nrows=1, ncols=1, figsize=(8*cm, 1*cm))
ax.bar(df_mAFiA_site_level['chromStart'], df_mAFiA_site_level['modRatio'])
ax.axhline(y=50, linestyle='--', c='r', alpha=0.1)
ax.set_xlim([x_start, x_end])
ax.set_ylim([0, 100])
ax.set_xticks(x_ticks)
ax.set_xlabel('Position (nts)')
ax.set_ylabel('Mod. Ratio')
fig_site_level.savefig(os.path.join(img_out, f'site_level_mAFiA.{FMT}'), **fig_kwargs)

with open(os.path.join(source_data_dir, 'source_data_Figure_1h.tsv'), 'w') as fout:
    fout.write('Figure 1h\n')
    fout.write('\n\t' + 'Position' + '\t')
    fout.write('\t'.join([str(x) for x in df_mAFiA_site_level['chromStart']]) + '\n')

    fout.write('\t' + 'Mod Ratio' + '\t')
    fout.write('\t'.join([str(x) for x in df_mAFiA_site_level['modRatio']]) + '\n')

########################################################################################################################
### vlnplot ############################################################################################################
########################################################################################################################
num_rows = 1
num_cols = 1
xtick_pos = [1, 2]
colors = ['pink', 'cyan']
cycles = ['even', 'odd']
fig_vlnplot, ax = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(5*cm, 5*cm))
mod_probs = [mAFiA_even_odd_probs[cycle] for cycle in cycles]
parts = ax.violinplot(mod_probs, xtick_pos, showmeans=False, showmedians=False, showextrema=False)
for pc, color in zip(parts['bodies'], colors):
    pc.set_facecolor(color)
    pc.set_edgecolor('black')
    pc.set_alpha(1)

    ax.set_xticks(xtick_pos, cycles)

    # if row_ind==(num_rows-1):
        # ax.set_xlabel('Aligned pattern')
    # else:
    #     ax.set_xticks([])
    ax.set_yticks(np.arange(0, 1.01, 0.5))
    # ax.set_ylabel(ds_names[test_ds])

    ax.set_xlim([0.5, 2.5])
    ax.set_ylim([-0.05, 1.05])
fig_vlnplot.savefig(os.path.join(img_out, f'vlnplot_mAFiA.{FMT}'), **fig_kwargs)

with open(os.path.join(source_data_dir, 'source_data_Figure_1i.tsv'), 'w') as fout:
    fout.write('Figure 1i\n\n')
    fout.write('\t' + 'P(m6A)' + '\t')
    fout.write('\t'.join([str(x * 0.01) for x in range(100)]) + '\n')
    for ind, this_cycle in enumerate(cycles):
        fout.write('\t' + this_cycle + '\t')
        out_data = np.histogram(mod_probs[ind], bins=100)[0]
        fout.write('\t'.join([str(x) for x in out_data]) + '\n')

########################################################################################################################
### precision-recall curve #############################################################################################
########################################################################################################################
def calc_prc(pattern_modProbs):
    thresholds = np.arange(0, 1.0, 0.01)
    prec = []
    rec = []
    for thresh in thresholds:
        tp = (pattern_modProbs['even'] > thresh).sum()
        fa = (pattern_modProbs['odd'] > thresh).sum()
        all_pos = len(pattern_modProbs['even'])

        if (tp==0) and (fa==0):
            break

        prec.append(tp / (tp + fa))
        rec.append(tp / all_pos)
    prec.append(1.0)
    rec.append(0.0)
    area = auc(rec[::-1], prec[::-1])
    return prec, rec, area

mAFiA_precisions, mAFiA_recalls, mAFiA_auprc = calc_prc(mAFiA_even_odd_probs)
CHEUI_precisions, CHEUI_recalls, CHEUI_auprc = calc_prc(CHEUI_even_odd_probs)
m6Anet_precisions, m6Anet_recalls, m6Anet_auprc = calc_prc(m6Anet_even_odd_probs)

fig_prc, ax = plt.subplots(nrows=1, ncols=1, figsize=(5*cm, 5*cm))
ax.plot(mAFiA_recalls, mAFiA_precisions, label=f'mAFiA, {mAFiA_auprc:.2f}')
ax.plot(CHEUI_recalls, CHEUI_precisions, label=f'CHEUI, {CHEUI_auprc:.2f}')
ax.plot(m6Anet_recalls, m6Anet_precisions, label=f'm6Anet, {m6Anet_auprc:.2f}')
ax.set_xticks(np.arange(0, 1.01, 0.5))
ax.set_yticks(np.arange(0, 1.01, 0.5))
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.legend(loc='lower left')

fig_prc.savefig(os.path.join(img_out, f'prc_cf_CHEUI_m6Anet.{FMT}'), **fig_kwargs)

method_recall_precision = {}
method_recall_precision['mAFiA'] = {}
method_recall_precision['mAFiA']['recall'] = mAFiA_recalls
method_recall_precision['mAFiA']['precision'] = mAFiA_precisions
method_recall_precision['CHEUI'] = {}
method_recall_precision['CHEUI']['recall'] = CHEUI_recalls
method_recall_precision['CHEUI']['precision'] = CHEUI_precisions
method_recall_precision['m6Anet'] = {}
method_recall_precision['m6Anet']['recall'] = m6Anet_recalls
method_recall_precision['m6Anet']['precision'] = m6Anet_precisions
with open(os.path.join(source_data_dir, 'source_data_Figure_1j.tsv'), 'w') as fout:
    fout.write('Figure 1j\n')
    for this_method in ['mAFiA', 'CHEUI', 'm6Anet']:
        fout.write('\n\t' + this_method + '\n')
        for this_label in ['recall', 'precision']:
            fout.write('\t' + this_label + '\t')
            fout.write('\t'.join([str(x) for x in method_recall_precision[this_method][this_label][::-1]]) + '\n')

########################################################################################################################
### construct probality mat ############################################################################################
########################################################################################################################
# min_pos = 6
# extent = 80
# min_cycles = 2
#
# for vmin in np.arange(0.1, 1, 0.1):
#
#     ds_mat_modProb = {}
#     for test_ds in test_datasets:
#         refPos_modProbs = ds_refPos_modProbs[test_ds]
#         # mod_probs = [pair for pair in ds_pattern_modProbs[test_ds][c_pattern] if np.max(pair)>=vmin]
#         # print(test_ds, c_pattern, len(mod_probs), np.mean(np.vstack(mod_probs), axis=0))
#
#         sample_refPos_modProbs = list(refPos_modProbs.values())
#         mat_mod_prob = np.zeros([len(refPos_modProbs), extent])
#         i = 0
#         for pos_probs in sample_refPos_modProbs:
#             refPos = np.array([tup[0] for tup in pos_probs])
#             modProbs = np.array([tup[1] for tup in pos_probs])
#             if np.sum(modProbs>vmin)<min_cycles:
#                 continue
#
#             first_high_ind = np.where(modProbs>vmin)[0][0]
#             # if first_high_ind==0:
#             #     filtered_refPos = refPos
#             #     filtered_modProbs = modProbs
#             #     shifted_refPos = filtered_refPos - np.min(filtered_refPos) + min_pos + block_size
#             # else:
#             filtered_refPos = refPos[first_high_ind:]
#             filtered_modProbs = modProbs[first_high_ind:]
#             shifted_refPos = filtered_refPos - np.min(filtered_refPos) + min_pos
#
#             extent_mask = shifted_refPos<extent
#             shifted_refPos = shifted_refPos[extent_mask]
#             filtered_modProbs = filtered_modProbs[extent_mask]
#
#             # if np.max(shifted_pos)>=extent:
#             #     continue
#             # mat_mod_prob[i, filtered_pos] = filtered_probs
#             mat_mod_prob[i, shifted_refPos] = filtered_modProbs
#             i+=1
#
#         mat_mod_prob = mat_mod_prob[~np.all(mat_mod_prob==0, axis=1)]
#         ds_mat_modProb[test_ds] = mat_mod_prob
#         print(f"{test_ds}: {len(mat_mod_prob)} samples")
#
#     ########################################################################################################################
#     ### visualize ##########################################################################################################
#     ########################################################################################################################
#     num_samples = 30
#
#     # num_samples = 10
#     # ytick_pos = [0, 4, 9]
#
#     fig_width = 5*cm
#     fig_height = 5*cm
#
#     num_rows = 2
#     num_cols = 1
#
#     # cax_yticks = np.linspace(vmin, 1, 3)
#     xticks = np.arange(min_pos, extent, block_size)
#     ytick_pos = np.concatenate([[0], np.arange(num_samples//2-1, num_samples, num_samples//2)])
#     yticks = [f'read {i+1}' for i in ytick_pos]
#
#     fig_single_read = plt.figure(figsize=(fig_width, fig_height))
#     for subplot_ind, test_ds in enumerate(test_datasets):
#         full_mat = ds_mat_modProb[test_ds]
#         # sel_rows = np.arange(num_samples)
#         sel_rows = np.argpartition(np.max(np.abs(np.diff(full_mat[:, np.arange(min_pos, extent, block_size)], axis=1)), axis=1), -num_samples)[-num_samples:]
#         display_mat = full_mat[sel_rows]
#         ax = fig_single_read.add_subplot(num_rows, num_cols, subplot_ind+1)
#         im = ax.imshow(display_mat, vmin=vmin, vmax=1, cmap='plasma')
#         # ax.set_xticks(xticks)
#         if subplot_ind==(num_rows-1):
#             ax.set_xticks(xticks)
#             ax.set_xlabel('Aligned Position')
#         else:
#             ax.set_xticks([])
#         ax.set_yticks([])
#         ax.set_yticks(ytick_pos, yticks)
#         ax.set_xlim([-0.5, extent-0.5])
#         # ax.set_title(ds_names[test_ds])
#
#     fig_single_read.tight_layout(rect=[0.1, 0.1, 0.8, 0.9])
#     fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
#     cax = fig_single_read.add_axes([0.85, 0.3, 0.03, 0.4])
#     plt.colorbar(im, cax=cax)
#     # plt.yticks([0.8, 0.9, 1.0])
#     plt.yticks(np.linspace(vmin, 1.0, 3))
#     plt.ylabel('P(m6A)', rotation=-90, labelpad=10)
#
#     fig_single_read.savefig(os.path.join(img_out, f'single_read_vmin{vmin:.1f}.{FMT}'), **fig_kwargs)
#     plt.close('all')
#
#     ########################################################################################################################
#     ### bar plots ##########################################################################################################
#     ########################################################################################################################
#     thresh_modProb = vmin
#     fig_barplot = plt.figure(figsize=(fig_width, fig_height))
#     ds_contrast = {}
#     for subplot_ind, test_ds in enumerate(test_datasets):
#         x_pos = np.arange(min_pos, extent, block_size)
#         prob_mat = ds_mat_modProb[test_ds][:, x_pos]
#         avg_mod_ratio = np.sum(prob_mat>=thresh_modProb, axis=0) / np.sum(prob_mat>0, axis=0)
#         ds_contrast[test_ds] = np.mean(avg_mod_ratio[[2, 4]] - avg_mod_ratio[[1, 3]])
#
#         ax = fig_barplot.add_subplot(num_rows, num_cols, subplot_ind+1)
#         ax.bar(x_pos, avg_mod_ratio, width=5, label=f'{ds_names[test_ds]}')
#         ax.axhline(y=0.5, c='r', linestyle='--', alpha=0.5)
#         if subplot_ind==(num_rows-1):
#             ax.set_xticks(xticks)
#             ax.set_xlabel('Aligned Position')
#         else:
#             ax.set_xticks([])
#         ax.set_ylim([0.0, 1.05])
#         ax.set_ylabel(f'Mod. Ratio')
#         # if subplot_ind==1:
#         #     ax.set_ylabel(f'Mod. Ratio')
#             # ax.set_ylabel(f'% NTs with P(m6A)>={thresh_modProb:.1f}')
#             # ax.set_ylabel('Mean P(m6A)')
#         ax.legend(loc='upper right')
#         # ax.set_yticks([])
#         # ax.set_yticks(ytick_pos, yticks)
#         # ax.set_title(ds_names[test_ds])
#     full_contrast = ds_contrast['WUE_batch3_AB_BA'] - ds_contrast['WUE_batch3_ABBA']
#     fig_barplot.savefig(os.path.join(img_out, f'barplot_vmin{vmin:.1f}_contrast{full_contrast:.2f}.{FMT}'), **fig_kwargs)
#     plt.close('all')
#
#     ########################################################################################################################
#     ### block distances ####################################################################################################
#     ########################################################################################################################
#     ds_diff = {}
#     fig_dist_freq = plt.figure(figsize=(fig_width, fig_height))
#     for subplot_ind, test_ds in enumerate(test_datasets):
#         mat_thresh = (ds_mat_modProb[test_ds]>=thresh_modProb)
#         i_ind, j_ind = np.where(mat_thresh)
#         dists = []
#         for i in np.unique(i_ind):
#             js = j_ind[i_ind == i]
#             dists.extend(list(np.diff(js)))
#         dist_counts = Counter(dists).most_common()
#         # block_dists = [int(tup[0]/block_size) for tup in dist_counts]
#         nt_dists = [int(tup[0]) for tup in dist_counts]
#         counts = [tup[1] for tup in dist_counts]
#         sorted_counts = np.array(counts)[np.argsort(nt_dists)]
#         freq = sorted_counts / np.sum(sorted_counts)
#         # block_dists = np.array(block_dists)[np.argsort(block_dists)]
#         nt_dists = np.array(nt_dists)[np.argsort(nt_dists)]
#
#         ds_diff[test_ds] = np.sum((freq[np.isin(nt_dists, [26, 52])] - freq[np.isin(nt_dists, [13, 39])]))
#
#         ax = fig_dist_freq.add_subplot(num_rows, num_cols, subplot_ind+1)
#         ax.bar(np.arange(len(nt_dists)), freq, width=0.5, label=ds_names[test_ds])
#         if subplot_ind==(num_rows-1):
#             ax.set_xticks(np.arange(len(nt_dists)), nt_dists)
#             ax.set_xlabel('m6A Distance')
#         else:
#             ax.set_xticks([])
#         ax.set_ylabel('Norm. Frequency')
#         ax.set_xlim([-0.5, 3.5])
#         ax.set_ylim([0.0, 0.65])
#         ax.set_yticks(np.arange(0, 0.65, 0.2))
#         ax.legend(loc='upper right')
#         # ax.set_title(ds_names[test_ds])
#
#     inter_ds_diff = ds_diff['WUE_batch3_AB_BA'] - ds_diff['WUE_batch3_ABBA']
#     # high = ds_diff['WUE_batch3_AB_BA']
#     # low = ds_diff['WUE_batch3_ABBA']
#
#     fig_dist_freq.savefig(os.path.join(img_out, f'distFreq_vmin{vmin:.1f}_diff{inter_ds_diff:.2f}.{FMT}'), **fig_kwargs)
#     plt.close('all')