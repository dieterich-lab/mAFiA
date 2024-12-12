import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import os

thresh_conf = 80.0
img_out = '/home/adrian/img_out/Gmorah'
dataset = 'rep1'
gmorah_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/Gmorah/HepG2/{dataset}/chrALL.mAFiA.sites.bed'
ref_file = '/home/adrian/Data/NmMutSeq/HepG2_mRNA_WT_Gm_sites.bed'

df_gmorah = pd.read_csv(gmorah_file, sep='\t', dtype={'chrom': str})

df_ref = pd.read_csv(ref_file, sep='\t', dtype={'chrom': str})
df_ref.rename(columns={'score': 'modRatio'}, inplace=True)

df_merge = pd.merge(df_gmorah[df_gmorah['confidence']>=thresh_conf], df_ref, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=['_Gmorah', '_NmMutSeq'])

plt.figure(figsize=(5, 5))
plt.plot(df_merge['modRatio_NmMutSeq'], df_merge['modRatio_Gmorah'], 'o')
# plt.plot(df_merge['modRatio_HeLa'], df_merge['modRatio_Gmorah'], 'o')
# plt.plot(df_merge['modRatio_HeLa'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['modRatio_HepG2'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['frac_RiboMethSeq'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['frac_NmMutSeq'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['modRatio_NmMutSeq'], df_merge['modRatio'], 'o')
plt.xlim([-5, 105])
plt.ylim([-5, 105])
plt.xlabel('Nm-Mut-Seq', fontsize=12)
plt.ylabel('Gmorah', fontsize=12)
plt.title(f'HepG2 {dataset}\n{len(df_merge)} Gm Sites', fontsize=15)
plt.savefig(os.path.join(img_out, f'scatter_Gmorah_vs_NmMutSeq_HepG2_{dataset}_mRNA.png'), bbox_inches='tight')