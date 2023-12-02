import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

df_wt = pd.read_csv(
    '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH5_retrain/100_WT_0_IVT/chr6/mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str})
df_ivt = pd.read_csv(
    '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH5_retrain/Mettl3-KO/chr6/mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str})

# six_motifs = [
#     'GGACT',
#     'GGACA',
#     'GAACT',
#     'AGACT',
#     'GGACC',
#     'TGACT'
# ]

df_merged = pd.merge(df_wt, df_ivt, on=['chrom', 'chromStart', 'chromEnd', 'score', 'strand', 'ref5mer'], suffixes=['_wt', '_ivt'])
df_merged_sel = df_merged[~df_merged['ref5mer'].isin(['AAACA', 'TAACT', 'TAACC'])]

plt.figure(figsize=(5, 5))
plt.scatter(df_merged_sel['modRatio_wt'], df_merged_sel['modRatio_ivt'], alpha=0.5)
plt.xlim([-1, 101])
plt.ylim([-1, 101])
plt.xlabel('$S_{WT}$')
plt.ylabel('$S_{IVT}$')
plt.title('chr1 DRACH sites')

Counter(df_merged[(df_merged['modRatio_ivt']>50)]['ref5mer'])