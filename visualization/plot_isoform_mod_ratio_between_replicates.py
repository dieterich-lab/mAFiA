import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
# FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

rep1 = '100_WT_0_IVT'
rep2 = 'P2_WT'

res_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
df_rep1 = pd.read_csv(os.path.join(res_dir, f'isoform_modRatio_{rep1}.tsv'), sep='\t')
df_rep2 = pd.read_csv(os.path.join(res_dir, f'isoform_modRatio_{rep2}.tsv'), sep='\t')

img_out = os.path.join('/home/adrian/img_out/MAFIA')
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

cov = 100

df_merge = pd.merge(df_rep1, df_rep2, on=['Chr', 'Pos', 'Gene', 'Transcript ID'], suffixes=('_rep1', '_rep2'))
df_filtered = df_merge[df_merge['Num. Reads_rep2']>=cov]

corr = np.corrcoef(df_filtered['Mod. Ratio_rep1'], df_filtered['Mod. Ratio_rep2'])[0, 1]

ticks = np.arange(0, 1.01, 0.5)
fig_width = 5 * cm
fig_height = 5 * cm

plt.figure(figsize=(fig_width, fig_height))
plt.scatter(df_filtered['Mod. Ratio_rep1'], df_filtered['Mod. Ratio_rep2'], s=0.3, label=f'Corr. {corr:.2f}')
plt.xlim([0, 1.05])
plt.ylim([0, 1.05])
plt.xticks(ticks)
plt.yticks(ticks)
# plt.title(f'$\\geq${cov}X in rep2\n{len(df_filtered)} isoforms')
plt.legend(loc='upper left')
plt.xlabel('Isoform Mod. Ratio, Rep. 1')
plt.ylabel('Isoform Mod. Ratio, Rep. 2')
plt.savefig(os.path.join(img_out, f'isoform_modRatio_comparison_{rep1}_{rep2}_cov{cov}.{FMT}'), **fig_kwargs)