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
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
train_dataset = 'WUE_batches1+2'
test_datasets = [
    'M4_M5',
    'M4_M5star',
    'M4star_M5',
    'M4star_M5star',
]
ds_names = {
    'M4_M5' : '$M_{4}M_{5}$',
    'M4_M5star' : '$M_{4}M_{5}^{*}$',
    'M4star_M5' : '$M_{4}^{*}M_{5}$',
    'M4star_M5star' : '$M_{4}^{*}M_{5}^{*}$'
}
contig_patterns = [
    '44',
    '45',
    '54',
    '55'
]

out_pickle = os.path.join(results_dir, 'M4_M5_modProbs.pkl')

img_out = os.path.join(HOME, 'img_out/MAFIA', 'M4_M5_single_read_{}'.format(train_dataset))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

dfs = {}
for test_dataset in test_datasets:
    path = os.path.join(results_dir, 'res_train_{}_test_RL_{}.tsv'.format(train_dataset, test_dataset))
    dfs[test_dataset] = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})

### collect mod probs based on contig pattern ###
block_size = 21
block_center = 10
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
            contig_pattern = sub_df['contig'].values[0].split('_')[1]
            pattern_start_blocks = [iter.span()[0] for iter in re.finditer(c_pattern, contig_pattern)]
            for this_start_block in pattern_start_blocks:
                ref_pos_first = this_start_block * block_size + block_center
                ref_pos_second = ref_pos_first + block_size
                if np.isin([ref_pos_first, ref_pos_second], sub_df['q_pos'].values).all():
                    mod_prob_pairs.append((
                        sub_df[sub_df['q_pos']==ref_pos_first]['mod_prob'].values[0],
                        sub_df[sub_df['q_pos']==ref_pos_second]['mod_prob'].values[0],
                    ))
        print('Mean pair mod. probs:', np.mean(np.vstack(mod_prob_pairs), axis=0))
        ds_pattern_modProbs[test_ds][c_pattern] = mod_prob_pairs

### dump to pickle ###
with open(out_pickle, 'wb') as handle:
    pickle.dump(ds_pattern_modProbs, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open(out_pickle, 'rb') as handle:
#     ds_pattern_modProbs = pickle.load(handle)

### visualize patterns ###
num_samples = 40
min_pos = 10
extent = 40
vmin = 0.8
cax_yticks = np.arange(vmin, 1.01, 0.1)
num_rows = 4
num_cols = 4
fig_single_read = plt.figure(figsize=(10, 10))
for row_ind, test_ds in enumerate(test_datasets):
    for col_ind, c_pattern in enumerate(contig_patterns):
        mod_probs = ds_pattern_modProbs[test_ds][c_pattern]
        # print(test_ds, c_pattern, len(mod_probs), np.mean(np.vstack(mod_probs), axis=0))
        if num_samples>=len(mod_probs):
            sample_mod_probs = mod_probs
        else:
            sample_mod_probs = sample(mod_probs, num_samples)
        mat_mod_prob = np.zeros([num_samples, extent])
        for i, prob_pair in enumerate(sample_mod_probs):
            mat_mod_prob[i, min_pos] = prob_pair[0]
            mat_mod_prob[i, min_pos+block_size] = prob_pair[1]

            xtick_pos = [min_pos, min_pos+block_size]
            xticks = ['P{}'.format(p) for p in c_pattern]

            ax = fig_single_read.add_subplot(num_rows, num_cols, row_ind * num_cols + col_ind + 1)
            im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1)
            if row_ind==(num_rows-1):
                ax.set_xticks(xtick_pos, xticks, fontsize=8)
                ax.set_xlabel('Aligned pattern', fontsize=12)
            else:
                ax.set_xticks([])
            if col_ind==0:
                ax.set_ylabel(ds_names[test_ds], fontsize=12)
            ax.set_yticks([])
fig_single_read.tight_layout(rect=[0.1, 0.2, 0.8, 0.8])
fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
cax = fig_single_read.add_axes([0.75, 0.3, 0.05, 0.4])
plt.colorbar(im, cax=cax)
plt.ylabel('P(m6A)', fontsize=12, rotation=-90, labelpad=30)
plt.yticks(cax_yticks, cax_yticks)

fig_single_read.savefig(os.path.join(img_out, 'single_read_{}by{}_grid.png'.format(num_rows, num_cols)), bbox_inches='tight')