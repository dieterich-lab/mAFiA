import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'pdf'
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

col0_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Arabidopsis_thaliana/col0/chr1/mAFiA.sites.bed'
vir1_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Arabidopsis_thaliana/vir1/chr1/mAFiA.sites.bed'
# miclip_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/parker_miclip_sites.tds'

img_out = '/home/adrian/NCOMMS_revision/images/ARABIDOPSIS'
os.makedirs(img_out, exist_ok=True)

df_col0 = pd.read_csv(col0_file, sep='\t', dtype={'chrom': str})
df_vir1 = pd.read_csv(vir1_file, sep='\t', dtype={'chrom': str})

# df_miclip = pd.read_csv(miclip_file, sep='\t',
#                         usecols=[0, 1, 2, 10],
#                         names=['chrom', 'chromStart', 'chromEnd', 'score'],
#                         dtype={'chrom': str})
# df_merged = pd.merge(df_mAFiA, df_miclip, on=['chrom', 'chromStart', 'chromEnd'], suffixes=['_mAFiA', '_miClip'])

df_merged = pd.merge(df_col0, df_vir1, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_col0', '_vir1'])

fig_scatter = plt.figure(figsize=(3*cm, 3*cm))
plt.scatter(df_merged['modRatio_col0'], df_merged['modRatio_vir1'], s=0.5, alpha=0.5)
plt.xlim([-1, 101])
plt.ylim([-1, 101])
plt.xticks([0, 50, 100])
plt.yticks([0, 50, 100])
plt.xlabel('$S_{col0}$')
plt.ylabel('$S_{vir1}$')
fig_scatter.savefig(os.path.join(img_out, f'scatter_col0_vir1.{FMT}'), **fig_kwargs)