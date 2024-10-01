import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

transcriptome_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_abundance/transcript_abundance_TAC_merged_vs_SHAM_merged.tsv'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'

df_transcriptome = pd.read_csv(transcriptome_file, sep='\t')

plt.figure(figsize=(4*cm, 4*cm))
plt.hist(df_transcriptome['log2fc'], range=[-3, 3], bins=50)
plt.savefig(os.path.join(img_out, f'hist_log2fc_transcript_abundance.{FMT}'), **fig_kwargs)
