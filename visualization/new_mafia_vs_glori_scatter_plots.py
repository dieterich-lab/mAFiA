import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

df_res = pd.read_csv('/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/results/train_DRACH_test_100_WT_0_IVT/mAFiA.sites.bed.merged', sep='\t')
df_glori = pd.read_csv('/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
df_glori['chrom'] = [chr.lstrip('chr') for chr in df_glori['Chr']]
df_glori['chromStart'] = df_glori['Sites'] - 1
df_glori['modRatio'] = np.int32(np.round(df_glori['Ratio'] * 100.0))

df_merged = pd.merge(df_res, df_glori, on=['chrom', 'chromStart'], suffixes=('_pred', '_glori'))
df_merged_sel = df_merged[df_merged['Pvalue'] < 1E-10]
motif_counts = Counter(df_merged_sel['ref5mer']).most_common()
motifs = [pair[0] for pair in motif_counts]

plt.figure(figsize=(16, 8))
for ind, motif in enumerate(motifs):
    sub_df = df_merged_sel[df_merged_sel['ref5mer'] == motif]
    corr = np.corrcoef(sub_df['modRatio_glori'], sub_df['modRatio_pred'])[0, 1]
    plt.subplot(3, 6, ind + 1)
    plt.plot(sub_df['modRatio_glori'], sub_df['modRatio_pred'], '.')
    plt.title(f'{motif}, {corr:.2f}')
    plt.xlim([-5, 105])
    plt.ylim([-5, 105])
