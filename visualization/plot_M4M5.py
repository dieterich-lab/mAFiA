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
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
train_dataset = 'ISA-WUE'
test_datasets = [
    'ISA_run4_M4M5',
    'ISA_run4_M4M5star',
    'ISA_run4_M4starM5',
    'ISA_run4_M4starM5star',
]
ds_names = {
    'ISA_run4_M4M5' : '$M_{4}M_{5}$',
    'ISA_run4_M4M5star' : '$M_{4}M_{5}^{*}$',
    'ISA_run4_M4starM5' : '$M_{4}^{*}M_{5}$',
    'ISA_run4_M4starM5star' : '$M_{4}^{*}M_{5}^{*}$'
}
contig_patterns = [
    'M4S0-M4S0',
    'M4S0-M5S0',
    'M5S0-M4S0',
    'M5S0-M5S0'
]

thresh_q = 70
block_size = 21
block_center = 10

out_pickle = os.path.join(results_dir, 'ISA_run4_M4M5_modProbs_q{}.pkl'.format(thresh_q))

img_out = os.path.join(HOME, 'img_out/MAFIA', 'res_train_{}_test_ISA_run4_M4M5'.format(train_dataset))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### load mod. probs. or collect from scratch ###########################################################################
########################################################################################################################
if os.path.exists(out_pickle):
    with open(out_pickle, 'rb') as handle:
        ds_pattern_modProbs = pickle.load(handle)
else:
    dfs = {}
    for test_dataset in test_datasets:
        path = os.path.join(results_dir, 'res_train_{}_test_{}_q{}.tsv'.format(train_dataset, test_dataset, thresh_q))
        dfs[test_dataset] = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})

    ### collect mod probs based on contig pattern ###
    ds_pattern_modProbs = {}
    for test_ds in test_datasets:
        print('\nDataset {}'.format(test_ds))
        df = dfs[test_ds]
        unique_reads = df['read_id'].unique()
        ds_pattern_modProbs[test_ds] = {}
        for c_pattern in contig_patterns:
            mod_prob_pairs = []
            print('\nCollecting pattern {}...'.format(c_pattern))
            for this_read in tqdm(unique_reads):
                sub_df = df[df['read_id']==this_read]

                # contig_pattern = sub_df['contig'].values[0].split('_')[1]
                # pattern_start_blocks = [iter.span()[0] for iter in re.finditer(c_pattern, contig_pattern)]

                ### TODO ###
                contig_pattern = sub_df['contig'].values[0].lstrip('ISA-')
                pattern_start_blocks = [iter.span()[0] for iter in re.finditer(c_pattern, contig_pattern)]

                for this_start_block in pattern_start_blocks:
                    ref_pos_first = this_start_block * block_size + block_center
                    ref_pos_second = ref_pos_first + block_size
                    if np.isin([ref_pos_first, ref_pos_second], sub_df['ref_pos'].values).all():
                        mod_prob_pairs.append((
                            sub_df[sub_df['ref_pos']==ref_pos_first]['mod_prob'].values[0],
                            sub_df[sub_df['ref_pos']==ref_pos_second]['mod_prob'].values[0],
                        ))
            print('Mean pair mod. probs:', np.mean(np.vstack(mod_prob_pairs), axis=0))
            ds_pattern_modProbs[test_ds][c_pattern] = mod_prob_pairs

    ### dump to pickle ###
    with open(out_pickle, 'wb') as handle:
        pickle.dump(ds_pattern_modProbs, handle, protocol=pickle.HIGHEST_PROTOCOL)

########################################################################################################################
### visualize single-reads #############################################################################################
########################################################################################################################
sel_test_ds = ['ISA_run4_M4M5star', 'ISA_run4_M4starM5']
sel_patterns = ['45', '54']

num_samples = 40
min_pos = 10
extent = 40
vmin = 0.8
cax_yticks = np.arange(vmin, 1.01, 0.1)
num_rows = len(sel_test_ds)
num_cols = len(sel_patterns)
fig_single_read = plt.figure(figsize=(num_cols*4+2, num_rows*4))
for row_ind, test_ds in enumerate(sel_test_ds):
    for col_ind, c_pattern in enumerate(sel_patterns):
        mod_probs = ds_pattern_modProbs[test_ds][c_pattern]
        # print(test_ds, c_pattern, len(mod_probs), np.mean(np.vstack(mod_probs), axis=0))
        if num_samples>=len(mod_probs):
            sample_mod_probs = mod_probs
        else:
            # sample_mod_probs = sample(mod_probs, num_samples)
            sample_mod_probs = [mod_probs[i] for i in np.argpartition(np.max(np.vstack(mod_probs), axis=1), -num_samples)[-num_samples:]]
        mat_mod_prob = np.zeros([num_samples, extent])
        for i, prob_pair in enumerate(sample_mod_probs):
            mat_mod_prob[i, min_pos] = prob_pair[0]
            mat_mod_prob[i, min_pos+block_size] = prob_pair[1]

            xtick_pos = [min_pos, min_pos+block_size]
            xticks = ['P{}'.format(p) for p in c_pattern]

            ax = fig_single_read.add_subplot(num_rows, num_cols, row_ind * num_cols + col_ind + 1)
            im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1)
            if row_ind==(num_rows-1):
                ax.set_xticks(xtick_pos, xticks, fontsize=10)
                ax.set_xlabel('Aligned pattern', fontsize=15)
            else:
                ax.set_xticks([])
            if col_ind==0:
                ax.set_ylabel(ds_names[test_ds], fontsize=15)
            ax.set_yticks([])
fig_single_read.tight_layout(rect=[0.1, 0.2, 0.8, 0.8])
fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
cax = fig_single_read.add_axes([0.85, 0.2, 0.03, 0.6])
plt.colorbar(im, cax=cax)
plt.ylabel('P(m6A)', fontsize=15, rotation=-90, labelpad=30)
plt.yticks(cax_yticks, cax_yticks)

fig_single_read.savefig(os.path.join(img_out, 'single_read_grid_q{}.png'.format(thresh_q)), bbox_inches='tight')

########################################################################################################################
### violin plots #######################################################################################################
########################################################################################################################
xtick_pos = [1, 2]
fig_vlnplot, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(num_cols*4, num_rows*4))
for row_ind, test_ds in enumerate(sel_test_ds):
    for col_ind, c_pattern in enumerate(sel_patterns):
        mod_probs = np.vstack(ds_pattern_modProbs[test_ds][c_pattern])
        ax = axs[row_ind, col_ind]
        # ax.violinplot(mod_probs, xtick_pos, points=50, widths=0.8)
        ax.boxplot(mod_probs, showfliers=False, whis=0.5)
        ax.axhline(y=0.5, c='r', alpha=0.5)
        xticks = ['P{}'.format(p) for p in c_pattern]
        if row_ind==(num_rows-1):
            ax.set_xticks(xtick_pos, xticks, fontsize=10)
            ax.set_xlabel('Aligned pattern', fontsize=15)
        else:
            ax.set_xticks([])
        if col_ind==0:
            ax.set_yticks(np.arange(0, 1.01, 0.5), fontsize=10)
            ax.set_ylabel(ds_names[test_ds], fontsize=15)
        else:
            ax.set_yticks([])
        ax.set_xlim([0.5, 2.5])
        ax.set_ylim([-0.05, 1.05])

        # axs[row_ind, col_ind].set_title('Custom violinplot 1', fontsize=fs)
fig_vlnplot.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
fig_vlnplot.suptitle('P(m6A)', fontsize=25)
fig_vlnplot.savefig(os.path.join(img_out, 'vlnplot_q{}.png'.format(thresh_q)), bbox_inches='tight')
