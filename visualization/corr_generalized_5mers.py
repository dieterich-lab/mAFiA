import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

glori_file = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv'
res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/chrX_generalized_5mer/mAFiA.sites.bed'
img_out = '/home/adrian/img_out'

df_glori = pd.read_csv(glori_file, usecols=['Chr', 'Sites', 'Strand', 'Ratio', 'P_adjust'])
df_glori = df_glori[df_glori['Chr']=='chrX']
df_res = pd.read_csv(res_file, sep='\t')

train_motifs = list(df_res['train5mer'].unique())
train_motifs.sort(key=lambda x: (x[1], x[3], x[0], x[4]))
test_motifs = [
    m for m in list(df_res['ref5mer'].unique())
    if ((m[0] in ['A', 'G', 'T']) and (m[1] in ['A', 'G']) and (m[3]=='C') and (m[4] in ['A', 'C', 'T']))
]
test_motifs.sort(key=lambda x: (x[1], x[3], x[0], x[4]))
# test_motifs = train_motifs + [m for m in test_motifs if m not in train_motifs]

corr_mat = np.zeros((len(train_motifs), len(test_motifs)))

train_test_corr_num = []
for ind_train, this_train_motif in enumerate(train_motifs):
    for ind_test, this_test_motif in enumerate(test_motifs):
        sub_df_res = df_res[(df_res['ref5mer']==this_test_motif) & (df_res['train5mer']==this_train_motif)]
        sub_df_glori = df_glori.set_index('Sites').loc[sub_df_res['chromEnd']]
        glori_ratio = sub_df_glori['Ratio'].values
        P_adjust = sub_df_glori['P_adjust'].values

        corr = np.corrcoef(glori_ratio, sub_df_res['modRatio'])[0, 1]
        num_sites = len(sub_df_res)
        train_test_corr_num.append((this_train_motif, this_test_motif, corr, num_sites))

        corr_mat[ind_train, ind_test] = corr

num_test_sites = [[tup[3] for tup in train_test_corr_num if tup[1]==motif][0] for motif in test_motifs]
composite_test_motifs = [f'{x} ({y})' for (x, y) in zip(test_motifs, num_test_sites)]

### visualize ###
divider = 8
vmin = 0.8
vmax = 0.9
sub_mat = corr_mat[:, divider:]
plt.figure(figsize=(8, 5))
im = plt.imshow(sub_mat, cmap='plasma', vmin=vmin, vmax=vmax)
plt.xticks(np.arange(sub_mat.shape[1]), composite_test_motifs[divider:], rotation=-90)
plt.yticks(np.arange(sub_mat.shape[0]), train_motifs)
plt.xlabel('Test 5-mer', fontsize=12)
plt.ylabel('Train 5-mer', fontsize=12)
cb = plt.colorbar(im, fraction=0.04, pad=0.04)
cb.set_ticks(np.arange(vmin, vmax+0.01, 0.05))
cb.set_label('Correlation', fontsize=12, rotation=-90, labelpad=20)
plt.title('mAFiA-GLORI correlation by train / test 5-mers\nchr X', fontsize=15)
plt.savefig(os.path.join(img_out, 'corr_by_train_test_5mer_NGACN.png'), bbox_inches='tight')

vmin = 0.6
vmax = 0.8
sub_mat = corr_mat[:, :divider]
plt.figure(figsize=(8, 5))
im = plt.imshow(sub_mat, cmap='plasma', vmin=vmin, vmax=vmax)
plt.xticks(np.arange(sub_mat.shape[1]), composite_test_motifs[:divider], rotation=-90)
plt.yticks(np.arange(sub_mat.shape[0]), train_motifs)
plt.xlabel('Test 5-mer', fontsize=12)
plt.ylabel('Train 5-mer', fontsize=12)
cb = plt.colorbar(im, fraction=0.04, pad=0.04)
cb.set_ticks(np.arange(vmin, vmax+0.01, 0.05))
cb.set_label('Correlation', fontsize=12, rotation=-90, labelpad=20)
plt.title('mAFiA-GLORI correlation by train / test 5-mers\nchr X', fontsize=15)
plt.savefig(os.path.join(img_out, 'corr_by_train_test_5mer_NAACN.png'), bbox_inches='tight')
