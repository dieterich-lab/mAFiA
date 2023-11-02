import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

img_out = '/home/adrian/img_out/mAFiA'
os.makedirs(img_out, exist_ok=True)

df_res = pd.read_csv('/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/results/train_DRACH_test_100_WT_0_IVT/mAFiA.sites.bed.merged', sep='\t')
df_glori = pd.read_csv('/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
df_glori['chrom'] = [chr.lstrip('chr') for chr in df_glori['Chr']]
df_glori['chromStart'] = df_glori['Sites'] - 1
df_glori['modRatio'] = np.int32(np.round(df_glori['NormeRatio'] * 100.0))

df_merged = pd.merge(df_res, df_glori, on=['chrom', 'chromStart'], suffixes=('_pred', '_glori'))
df_merged_sel = df_merged[df_merged['P_adjust'] < 1E-10]
motif_counts = Counter(df_merged_sel['ref5mer']).most_common()
motifs = [pair[0] for pair in motif_counts]

total_num_sites = df_merged_sel.shape[0]
total_corr = np.corrcoef(df_merged_sel['modRatio_glori'], df_merged_sel['modRatio_pred'])[0, 1]

num_rows = 3
num_cols = 6

plt.figure(figsize=(20, 10))
for ind, motif in enumerate(motifs):
    sub_df = df_merged_sel[df_merged_sel['ref5mer'] == motif]
    corr = np.corrcoef(sub_df['modRatio_glori'], sub_df['modRatio_pred'])[0, 1]
    plt.subplot(num_rows, num_cols, ind + 1)
    plt.plot(sub_df['modRatio_glori'], sub_df['modRatio_pred'], '.')
    plt.title(f'{motif}, {corr:.2f}', x=0.26, y=1.0, pad=-15, backgroundcolor='black', color='white')
    plt.xlim([-5, 105])
    plt.ylim([-5, 105])

    plt.xticks(range(0, 101, 25))
    plt.yticks(range(0, 101, 25))

    if ind%num_cols == 0:
        plt.ylabel('mAFiA', fontsize=12)
    if ind>=(num_rows-1)*num_cols:
        plt.xlabel('GLORI', fontsize=12)
plt.suptitle(f'HEK293 WT\nTotal {total_num_sites} sites\nCorr. {total_corr:.2f}', fontsize=15)
plt.savefig(os.path.join(img_out, 'corr_mAFiA_GLORI_DRACH.png'), bbox_inches='tight')
plt.close('all')
