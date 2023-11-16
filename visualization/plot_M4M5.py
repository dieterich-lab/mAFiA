import os
HOME = os.path.expanduser('~')
import re
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import random
random.seed(0)
from random import sample
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics import auc

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 7
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/results'
train_dataset = 'ISA-WUE'
test_datasets = [
    # 'ISA_run4_M4M5',
    'ISA_run4_M4M5star',
    'ISA_run4_M4starM5',
    # 'ISA_run4_M4starM5star',
]
# ds_names = {
#     # 'ISA_run4_M4M5' : 'GGACC + TGACT',
#     'ISA_run4_M4M5star' : 'GGACC + TGm6ACT',
#     'ISA_run4_M4starM5' : 'GGm6ACC + TGACT',
#     # 'ISA_run4_M4starM5star' : 'GGm6ACC + TGm6ACT'
# }

ds_names = {
    # 'ISA_run4_M4M5' : 'GGACC + TGACT',
    'ISA_run4_M4M5star' : 'REP1',
    'ISA_run4_M4starM5' : 'REP2',
    # 'ISA_run4_M4starM5star' : 'GGm6ACC + TGm6ACT'
}

contig_patterns = [
    # 'M4S0-M4S0',
    'M4S0-M5S0',
    'M5S0-M4S0',
    # 'M5S0-M5S0'
]

# sequences = {
#     'M4S0' : 'CCCCAGCTGGACCGACTCAGA',
#     'M5S0' : 'GGGACCTCTGACTGCTCTGGG'
# }

contig_motifs = {
    'M4S0' : 'GGACC',
    'M5S0' : 'TGACT'
}

thresh_q = 70
block_size = 21
block_center = 10

img_out = os.path.join(HOME, 'img_out/NCOMMS_rev', 'res_train_{}_test_ISA_run4_M4M5_minimap'.format(train_dataset))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### load mod. probs. or collect from scratch ###########################################################################
########################################################################################################################
dfs = {}
for test_dataset in test_datasets:
    # path = os.path.join(results_dir, 'res_train_{}_test_{}_q{}.tsv'.format(train_dataset, test_dataset, thresh_q))
    path = os.path.join(results_dir, 'res_train_{}_test_{}_minimap.tsv'.format(train_dataset, test_dataset))
    this_df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})
    dfs[test_dataset] = this_df
    # dfs[test_dataset] = this_df[this_df['ref_motif']==this_df['pred_motif']]

### collect modProbs ###
ds_pattern_modProbs = {}
for test_ds in test_datasets:
    df = dfs[test_ds]
    ds_pattern_modProbs[test_ds] = {}
    for motif in list(contig_motifs.values()):
        ds_pattern_modProbs[test_ds][motif] = df[df['ref_motif']==motif]['mod_prob'].values

########################################################################################################################
### vlnplot ############################################################################################################
########################################################################################################################
num_rows = 2
num_cols = 1
xtick_pos = [1, 2]
colors = ['pink', 'cyan']
fig_vlnplot, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(3*cm, 5*cm))
for row_ind, test_ds in enumerate(test_datasets):
    mod_probs = [ds_pattern_modProbs[test_ds][motif] for motif in list(contig_motifs.values())]
    ax = axs[row_ind]
    parts = ax.violinplot(mod_probs, xtick_pos, showmeans=False, showmedians=False, showextrema=False)
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    ax.set_xticks(xtick_pos, contig_motifs.values())

    # if row_ind==(num_rows-1):
        # ax.set_xlabel('Aligned pattern')
    # else:
    #     ax.set_xticks([])
    ax.set_yticks(np.arange(0, 1.01, 0.5))
    # ax.set_ylabel(ds_names[test_ds])

    ax.set_xlim([0.5, 2.5])
    ax.set_ylim([-0.05, 1.05])
fig_vlnplot.savefig(os.path.join(img_out, f'vlnplot_q{thresh_q}.{FMT}'), **fig_kwargs)

########################################################################################################################
### parse CHEUI output #################################################################################################
########################################################################################################################
CHEUI_output_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/CHEUI/oligo'

CHEUI_pattern_modProbs = {}
for test_ds in test_datasets:
    df = pd.read_csv(os.path.join(CHEUI_output_dir, test_ds, 'm6A/read_level_m6A_predictions.txt'), sep='\t', names=['segment', 'mod_prob', 'rep'])
    CHEUI_pattern_modProbs[test_ds] = {}
    mod_probs = df['mod_prob'].values
    segs = np.array([seg.split('_')[3] for seg in df['segment'].values])
    CHEUI_pattern_modProbs[test_ds]['GGACC'] = mod_probs[segs=='CTGGACCGA']
    CHEUI_pattern_modProbs[test_ds]['TGACT'] = mod_probs[segs=='TCTGACTGC']

########################################################################################################################
### parse m6Anet output ################################################################################################
########################################################################################################################
m6Anet_output_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/m6Anet/oligo'

m6Anet_pattern_modProbs = {}
for test_ds in test_datasets:
    m6Anet_pattern_modProbs[test_ds] = {}
    df_indiv = pd.read_csv(os.path.join(m6Anet_output_dir, test_ds, 'data.indiv_proba.csv'))
    df_site = pd.read_csv(os.path.join(m6Anet_output_dir, test_ds, 'data.site_proba.csv'))
    df_indiv_central = df_indiv[df_indiv['transcript_position'] % 21 == 10]

    kmers = []
    for _, row in df_indiv_central.iterrows():
        sub_df = df_site[(df_site['transcript_id']==row['transcript_id']) * (df_site['transcript_position']==row['transcript_position'])]
        if len(sub_df)>1:
            print('Non-unique site!!!')
        kmers.append([str(val) for val in sub_df['kmer'].values])
    kmers = [xx for ll in kmers for xx in ll]
    df_indiv_central['kmer'] = kmers

    m6Anet_pattern_modProbs[test_ds]['GGACC'] = df_indiv_central[df_indiv_central['kmer']=='GGACC']['probability_modified'].values
    m6Anet_pattern_modProbs[test_ds]['TGACT'] = df_indiv_central[df_indiv_central['kmer']=='TGACT']['probability_modified'].values

########################################################################################################################
### precision-recall curve #############################################################################################
########################################################################################################################
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5

def calc_prc(pattern_modProbs):
    thresholds = np.arange(0, 1.0, 0.01)
    prec = []
    rec = []
    for thresh in thresholds:
        tp = (pattern_modProbs['ISA_run4_M4M5star']['TGACT'] > thresh).sum() + (pattern_modProbs['ISA_run4_M4starM5']['GGACC'] > thresh).sum()
        fa = (pattern_modProbs['ISA_run4_M4M5star']['GGACC'] > thresh).sum() + (pattern_modProbs['ISA_run4_M4starM5']['TGACT'] > thresh).sum()
        all_pos = len(pattern_modProbs['ISA_run4_M4M5star']['TGACT']) + len(pattern_modProbs['ISA_run4_M4starM5']['GGACC'])

        if (tp==0) and (fa==0):
            break

        prec.append(tp / (tp + fa))
        rec.append(tp / all_pos)
    prec.append(1.0)
    rec.append(0.0)
    area = auc(rec[::-1], prec[::-1])
    return prec, rec, area

precisions, recalls, auprc = calc_prc(ds_pattern_modProbs)
CHEUI_precisions, CHEUI_recalls, CHEUI_auprc = calc_prc(CHEUI_pattern_modProbs)
m6Anet_precisions, m6Anet_recalls, m6Anet_auprc = calc_prc(m6Anet_pattern_modProbs)

fig_prc, ax = plt.subplots(nrows=1, ncols=1, figsize=(5*cm, 5*cm))
ax.plot(recalls, precisions, label=f'mAFiA, {auprc:.2f}')
ax.plot(CHEUI_recalls, CHEUI_precisions, label=f'CHEUI, {CHEUI_auprc:.2f}')
ax.plot(m6Anet_recalls, m6Anet_precisions, label=f'm6Anet, {m6Anet_auprc:.2f}')
ax.set_xticks(np.arange(0, 1.01, 0.5))
ax.set_yticks(np.arange(0, 1.01, 0.5))
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.legend(loc='lower left')

fig_prc.savefig(os.path.join(img_out, f'prc_cf_CHEUI_m6Anet.{FMT}'), **fig_kwargs)


########################################################################################################################
### collect mod probs based on contig pattern ###
# ds_pattern_modProbs = {}
# for test_ds in test_datasets:
#     print('\nDataset {}'.format(test_ds))
#     df = dfs[test_ds]
#     unique_reads = df['read_id'].unique()
#     ds_pattern_modProbs[test_ds] = {}
#     for this_read in tqdm(unique_reads):
#         sub_df = df[df['read_id']==this_read]
#         # contig_pattern = sub_df['contig'].values[0].split('_')[1]
#         contig_pattern = sub_df['contig'].values[0].lstrip('ISA-')
#
#         for c_pattern in contig_patterns:
#             # print('\nCollecting pattern {}...'.format(c_pattern))
#             if c_pattern not in ds_pattern_modProbs[test_ds].keys():
#                 ds_pattern_modProbs[test_ds][c_pattern] = []
#
#             # pattern_start_blocks = [it.span()[0] for it in re.finditer(c_pattern, contig_pattern)]
#             pattern_start_blocks = [it.span()[0] for it in re.finditer(c_pattern, contig_pattern)]
#
#             mod_prob_pairs = []
#             for this_start_block in pattern_start_blocks:
#                 ref_pos_first = this_start_block * block_size + block_center
#                 ref_pos_second = ref_pos_first + block_size
#                 if (np.isin([ref_pos_first, ref_pos_second], sub_df['ref_pos'].values).all())\
#                         and (np.abs((sub_df[sub_df['ref_pos']==ref_pos_second]['read_pos'].values-sub_df[sub_df['ref_pos']==ref_pos_first]['read_pos'].values)[0]-block_size)<=1):
#                     mod_prob_pairs.append((
#                         sub_df[sub_df['ref_pos']==ref_pos_first]['mod_prob'].values[0],
#                         sub_df[sub_df['ref_pos']==ref_pos_second]['mod_prob'].values[0],
#                     ))
#             ds_pattern_modProbs[test_ds][c_pattern].extend(mod_prob_pairs)

### dump to pickle ###
# with open(out_pickle, 'wb') as handle:
#     pickle.dump(ds_pattern_modProbs, handle, protocol=pickle.HIGHEST_PROTOCOL)

########################################################################################################################
### visualize single-reads #############################################################################################
########################################################################################################################
# sel_test_ds = ['ISA_run4_M4M5star', 'ISA_run4_M4starM5']
# sel_patterns = ['M4S0-M5S0', 'M5S0-M4S0']
#
# fig_width = (5 + 1)*cm
# fig_height = 5*cm
#
# num_samples = 40
# ytick_pos = [0, 19, 39]
# yticks = [f'read {i+1}' for i in ytick_pos]
# min_pos = 10
# extent = 2 * block_size
# vmin = 0.8
# cax_yticks = np.linspace(vmin, 1, 3)
# num_rows = len(sel_test_ds)
# num_cols = len(sel_patterns)
#
# fig_single_read = plt.figure(figsize=(fig_width, fig_height))
# for row_ind, test_ds in enumerate(sel_test_ds):
#     for col_ind, c_pattern in enumerate(sel_patterns):
#         # mod_probs = ds_pattern_modProbs[test_ds][c_pattern]
#         # mod_probs = [pair for pair in ds_pattern_modProbs[test_ds][c_pattern] if np.max(pair)>=vmin]
#         mod_probs = [pair for pair in ds_pattern_modProbs[test_ds][c_pattern] if ((np.max(pair)-np.min(pair))>=0.6)]
#
#         # print(test_ds, c_pattern, len(mod_probs), np.mean(np.vstack(mod_probs), axis=0))
#         xtick_pos = [min_pos, min_pos+block_size]
#         # xticks = ['P{}'.format(p) for p in c_pattern]
#         # xticks = ''.join([sequences[contig] for contig in c_pattern.split('-')])
#         xticks = [contig_motifs[contig] for contig in c_pattern.split('-')]
#         if num_samples>=len(mod_probs):
#             sample_mod_probs = mod_probs
#         else:
#             sample_mod_probs = sample(mod_probs, num_samples)
#         mat_mod_prob = np.zeros([num_samples, extent])
#         for i, prob_pair in enumerate(sample_mod_probs):
#             mat_mod_prob[i, min_pos] = prob_pair[0]
#             mat_mod_prob[i, min_pos+block_size] = prob_pair[1]
#
#         ax = fig_single_read.add_subplot(num_rows, num_cols, row_ind * num_cols + col_ind + 1)
#         im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1, cmap='plasma')
#         if row_ind==(num_rows-1):
#             ax.set_xticks(xtick_pos, xticks)
#             # ax.set_xlabel('Aligned pattern', fontsize=15)
#         else:
#             ax.set_xticks([])
#         if col_ind==0:
#             ax.set_yticks(ytick_pos, yticks)
#             ax.set_ylabel(ds_names[test_ds])
#         else:
#             ax.set_yticks([])
# # fig_single_read.tight_layout(rect=[0.1, 0.1, 0.8, 0.9])
# fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
# cax = fig_single_read.add_axes([0.85, 0.3, 0.03, 0.4])
# plt.colorbar(im, cax=cax)
# # plt.ylabel('P(m6A)', fontsize=15, rotation=-90, labelpad=30)
# plt.yticks(cax_yticks, cax_yticks)
#
# fig_single_read.savefig(os.path.join(img_out, f'single_read_grid_q{thresh_q}.{FMT}'), **fig_kwargs)

########################################################################################################################
### box plots ##########################################################################################################
########################################################################################################################
# xtick_pos = [1, 2]
# fig_boxplot, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(5*cm, 5*cm))
# for row_ind, test_ds in enumerate(sel_test_ds):
#     for col_ind, c_pattern in enumerate(sel_patterns):
#         mod_probs = np.vstack(ds_pattern_modProbs[test_ds][c_pattern])
#         ax = axs[row_ind, col_ind]
#         # ax.violinplot(mod_probs, xtick_pos, points=50, widths=0.8)
#         ax.boxplot(mod_probs, showfliers=False, whis=0.5)
#         # ax.axhline(y=0.5, c='r', alpha=0.5)
#         # xticks = ['P{}'.format(p) for p in c_pattern]
#         xticks = [contig_motifs[contig] for contig in c_pattern.split('-')]
#
#         if row_ind==(num_rows-1):
#             ax.set_xticks(xtick_pos, xticks)
#             # ax.set_xlabel('Aligned pattern')
#         else:
#             ax.set_xticks([])
#         if col_ind==0:
#             ax.set_yticks(np.arange(0, 1.01, 0.5))
#             ax.set_ylabel(ds_names[test_ds])
#         else:
#             ax.set_yticks([])
#         ax.set_xlim([0.5, 2.5])
#         ax.set_ylim([-0.05, 1.05])
#
#         # axs[row_ind, col_ind].set_title('Custom violinplot 1', fontsize=fs)
# fig_boxplot.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
# # fig_boxplot.suptitle('P(m6A)')
# fig_boxplot.savefig(os.path.join(img_out, f'boxplot_q{thresh_q}.{FMT}'), **fig_kwargs)
