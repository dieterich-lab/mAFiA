import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import os

img_out = '/home/adrian/img_out/Gmorah'

# gmorah_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HeLa_rRNA/Gmorah/mAFiA.sites.bed'
# gmorah_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HEK293/Gmorah_HeLa_HepG2_merge/mAFiA.sites.bed'
# gmorah_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HEK293/Gmorah_HeLa/mAFiA.sites.bed'
gmorah_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HEK293_P2/Gmorah/mAFiA.sites.bed'
df_gmorah = pd.read_csv(gmorah_file, sep='\t')

# ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Nm-Mut-seq_Supp_Tables.xlsx'
# df_ref = pd.read_excel(ref_file, sheet_name='S1_HeLa_known sites_rRNA_WT', skiprows=[0, 1],
#                        usecols=['chr', 'position', 'motif', 'Fraction (%)', 'ref_Fraction% (SILNAS)', 'ref_Fraction% (RiboMethseq)'])
# df_ref = df_ref.rename(columns={
#     'chr': 'chrom',
#     'position': 'origPosition',
#     'motif': 'ref5mer',
#     'Fraction (%)': 'modRatio_NmMutSeq',
#     'ref_Fraction% (SILNAS)': 'modRatio_SILNAS',
#     'ref_Fraction% (RiboMethseq)': 'modRatio_RiboMethSeq'
# })

# ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/Gm_sites_mRNA_HeLa.bed'
# df_ref = pd.read_csv(ref_file, sep='\t')

# df_merge = pd.merge(df_gmorah, df_ref, on=['chrom', 'origPosition', 'ref5mer'], suffixes=['_gmorah', '_ref'])
# df_merge = pd.merge(df_gmorah, df_ref, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_Gmorah', '_HeLa'])
df_merge = df_gmorah

plt.figure(figsize=(5, 5))
# plt.plot(df_merge['modRatio_HeLa'], df_merge['modRatio_Gmorah'], 'o')
plt.plot(df_merge['modRatio_HeLa'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['modRatio_HepG2'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['frac_RiboMethSeq'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['frac_NmMutSeq'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['modRatio_NmMutSeq'], df_merge['modRatio'], 'o')
plt.xlim([-5, 105])
plt.ylim([-5, 105])
# plt.xlabel('RiboMethSeq', fontsize=12)
# plt.xlabel('NmMutSeq')
# plt.xlabel('HeLa NmMutSeq', fontsize=12)
plt.xlabel('HepG2 NmMutSeq', fontsize=12)
plt.ylabel('HEK293 Gmorah', fontsize=12)
# plt.title('HeLa rRNA 18S 28S Gm', fontsize=15)
plt.title('HeLa-HepG2 Gm Sites', fontsize=15)
# plt.savefig(os.path.join(img_out, 'scatter_Gmorah_vs_NmMutSeq_rRNA.png'), bbox_inches='tight')
plt.savefig(os.path.join(img_out, 'scatter_Gmorah_HEK293_P2_vs_NmMutSeq_HepG2_mRNA.png'), bbox_inches='tight')