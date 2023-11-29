import pandas as pd
import matplotlib.pyplot as plt

wt_file = '/home/adrian/NCOMMS_revision/source_data/YEAST/WT/chrI/mAFiA.sites.bed'
ko_file = '/home/adrian/NCOMMS_revision/source_data/YEAST/IME4_KO/chrI/mAFiA.sites.bed'

df_wt = pd.read_csv(wt_file, sep='\t')
df_ko = pd.read_csv(ko_file, sep='\t')
df_merged = pd.merge(df_wt, df_ko, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_wt', '_ko'])


plt.figure(figsize=(5, 5))
plt.plot(df_merged['modRatio_wt'], df_merged['modRatio_ko'], 'bo', alpha=0.5)