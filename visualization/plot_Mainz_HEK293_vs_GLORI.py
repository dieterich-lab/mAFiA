import pandas as pd
import matplotlib.pyplot as plt

glori = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv'
dmso = '/home/adrian/Data/TRR319_RMaP/Project_B01/Adrian/JK_HEK293_DMSO_merged/mAFiA/mAFiA.sites.bed'
stm = '/home/adrian/Data/TRR319_RMaP/Project_B01/Adrian/JK_HEK293_STM2457_merged/mAFiA/mAFiA.sites.bed'

min_cov = 50

dict_chr = {f'chr{i}': i for i in range(1, 25)}

df_glori = pd.read_csv(glori, usecols=['Chr', 'Sites', 'Strand', 'Gene', 'Transcript', 'Ratio', 'P_adjust'])
df_dmso = pd.read_csv(dmso, sep='\t')
df_stm = pd.read_csv(stm, sep='\t')

df_glori['chrom'] = df_glori['Chr'].apply(lambda x: x.lstrip('chr'))
df_glori = df_glori.rename(columns={'Sites' : 'chromEnd'})
df_glori['modRatio_GLORI'] = df_glori['Ratio'].apply(lambda x: int(x*100))

df_mafia = pd.merge(df_dmso, df_stm, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=['_DMSO', '_STM'])
df_glori_mafia = pd.merge(df_glori, df_mafia, on=['chrom', 'chromEnd'])
df_glori_mafia_filtered = df_glori_mafia[(df_glori_mafia['coverage_DMSO']>=min_cov) * (df_glori_mafia['coverage_STM']>=min_cov)]

plt.figure(figsize=(5, 5))
plt.plot(df_glori_mafia_filtered['modRatio_GLORI'], df_glori_mafia_filtered['modRatio_DMSO'], 'r.', label='DMSO')
plt.plot(df_glori_mafia_filtered['modRatio_GLORI'], df_glori_mafia_filtered['modRatio_STM'], 'b.', label='STM')
plt.xlabel('GLORI')
plt.ylabel(f'mAFiA, cov>={min_cov}')
plt.xlim([-1, 101])
plt.ylim([-1, 101])
plt.title('Mainz HEK293')
plt.legend(loc='best')
plt.savefig('/home/adrian/Data/TRR319_RMaP/Project_B01/Adrian/DMSO_STM_GLORI.png', bbox_inches='tight')
plt.close('all')