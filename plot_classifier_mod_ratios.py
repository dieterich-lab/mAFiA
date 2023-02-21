import os
import pandas as pd
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

img_out = '/home/adrian/img_out'

tsv_dir = '/home/adrian/inference/rRNA/METTL5_mixing_viterbi/'
mod = 'm6A'
stop = 1832
mixes = [
    'IVT',
    '0WT_100KO',
    '50WT_50KO',
    '100WT_0KO',
    'WT'
]

# tsv_dir = '/home/adrian/inference/rRNA/METTL5_mixing'
# mod = 'm6A'
# stop = 1832
# mixes = [
#     '0WT_100KO',
#     '25WT_75KO',
#     '50WT_50KO',
#     '75WT_25KO',
#     '100WT_0KO',
# ]

classifier = 'log_reg'

mix_mod_ratio = {}
mix_num_features = {}
for this_mix in mixes:
    df = pd.read_csv(os.path.join(tsv_dir, '{}_mod_ratio_{}.tsv'.format(classifier, this_mix)), sep='\t')
    # df = pd.read_csv(os.path.join(tsv_dir, 'svm_opt_thresh_mod_ratio_{}.tsv'.format(this_mix)), sep='\t')
    df_mod = df[(df['stop'] == stop) & (df['mod'] == mod)]
    mod_ratio = df_mod['mod_ratio'].values[0]
    num_test_features = df_mod['num_test_features'].values[0]
    mix_mod_ratio[this_mix] = mod_ratio
    mix_num_features[this_mix] = num_test_features

plt.figure(figsize=(8, 8))
plt.plot(range(len(mixes)), mix_mod_ratio.values(), 'o', mfc='None')
plt.xlim([-0.5, len(mixes)-0.5])
plt.xticks(range(len(mixes)), mixes)
plt.ylim([-10, 110])
plt.xlabel('Mix', fontsize=15)
plt.ylabel('Pred. mod. ratio', fontsize=15)
plt.title('{} {}'.format(mod, stop), fontsize=20)

plt.savefig(os.path.join(img_out, '{}_mod_ratios_{}_{}.png'.format(classifier, mod, stop)), bbox_inches='tight')
plt.close('all')