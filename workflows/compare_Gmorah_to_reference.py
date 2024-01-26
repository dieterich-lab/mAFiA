import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import os

img_out = '/home/adrian/img_out/Gmorah'

gmorah_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HeLa_rRNA/Gmorah/mAFiA.sites.bed'
df_gmorah = pd.read_csv(gmorah_file, sep='\t')

ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Nm-Mut-seq_Supp_Tables.xlsx'
df_ref = pd.read_excel(ref_file, sheet_name='S1_HeLa_known sites_rRNA_WT', skiprows=[0, 1],
                       usecols=['chr', 'position', 'motif', 'Fraction (%)', 'ref_Fraction% (SILNAS)', 'ref_Fraction% (RiboMethseq)'])
df_ref = df_ref.rename(columns={
    'chr': 'chrom',
    'position': 'origPosition',
    'motif': 'ref5mer',
    'Fraction (%)': 'frac_NmMutSeq',
    'ref_Fraction% (SILNAS)': 'frac_SILNAS',
    'ref_Fraction% (RiboMethseq)': 'frac_RiboMethSeq'
})

df_merge = pd.merge(df_gmorah, df_ref, on=['chrom', 'origPosition', 'ref5mer'], suffixes=['_gmorah', '_ref'])

plt.figure(figsize=(5, 5))
plt.plot(df_merge['frac_RiboMethSeq'], df_merge['modRatio'], 'o')
# plt.plot(df_merge['frac_NmMutSeq'], df_merge['modRatio'], '.')
plt.xlim([-5, 105])
plt.ylim([-5, 105])
plt.xlabel('RiboMethSeq', fontsize=12)
# plt.xlabel('NmMutSeq')
plt.ylabel('Gmorah', fontsize=12)
plt.title('HeLa rRNA 18S 28S Gm', fontsize=15)
plt.savefig(os.path.join(img_out, 'scatter_Gmorah_vs_ref_HeLa_rRNA.png'), bbox_inches='tight')