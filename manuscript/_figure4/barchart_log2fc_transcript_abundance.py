import pandas as pd
import numpy as np
import os
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
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
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200, transparent=True)
#######################################################################

data_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_abundance'
transcriptome_file = os.path.join(data_dir, 'transcript_abundance_TAC_merged_vs_SHAM_merged.tsv')
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/_figure4'

df_transcriptome = pd.read_csv(transcriptome_file, sep='\t')

# plt.figure(figsize=(4*cm, 4*cm))
# plt.hist(df_transcriptome['log2fc'], range=[-3, 3], bins=50)
# plt.savefig(os.path.join(img_out, f'hist_log2fc_transcript_abundance.{FMT}'), **fig_kwargs)

num_all_genes = 2000
sel_genes = []
for this_gene in df_transcriptome.sort_values('norm_coverage_SHAM_merged', ascending=False)['gene']:
    if len(sel_genes) >= num_all_genes:
        break
    if this_gene[:3] != 'mt-':
        sel_genes.append(this_gene)
sel_df = df_transcriptome[df_transcriptome['gene'].isin(sel_genes)].sort_values('log2fc')

thresh_log2fc = 1
df_negative = sel_df[sel_df['log2fc'] < -thresh_log2fc]
df_positive = sel_df[sel_df['log2fc'] >= thresh_log2fc]

with open(os.path.join(data_dir,
                       f'up_down_regulated_transcripts_common{num_all_genes}_log2fc{thresh_log2fc}.txt'), 'w') as f_out:
    f_out.write('#up_regulated\n')
    for val in df_positive['gene'].values:
        f_out.write(f'{val}\n')
    f_out.write('#down_regulated\n')
    for val in df_negative['gene'].values:
        f_out.write(f'{val}\n')

xlim = [-3, 3]
sel_genes = list(df_negative['gene'].values) + list(df_positive['gene'].values)
num_genes = len(sel_genes)
plt.figure(figsize=(4*cm, 6*cm))
plt.barh(np.arange(len(df_negative)), df_negative['log2fc'], fc='r')
plt.barh(np.arange(len(df_negative), num_genes), df_positive['log2fc'], fc='b')
plt.gca().invert_yaxis()
plt.yticks(np.arange(num_genes), sel_genes)
# plt.xlim(xlim)
plt.xticks(np.linspace(*xlim, 5))
plt.axvline(x=0, c='gray', ls='--')
plt.ylim([-0.8, num_genes-1+0.8])
plt.savefig(
    os.path.join(img_out,
                 f'log2fc_transcript_abundance_{num_genes}_most_common_transcripts_log2fc{thresh_log2fc}.{FMT}'),
    **fig_kwargs
)