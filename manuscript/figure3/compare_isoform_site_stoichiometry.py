import os
import pandas as pd
from functools import reduce
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
#######################################################################
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import seaborn as sns
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
# sns.set_style('whitegrid')

thresh_coverage = 10

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

data_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/mouse_heart'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'

ds = 'HFpEF'
conditions = ['ctrl', 'HFpEF']
bambu_file = os.path.join(
    '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/DTU',
    'DTU_HFpEF_Laura_ENS102.tsv'
)

df_bambu = pd.read_csv(bambu_file, sep='\t')

gene = 'Arpp19'

transcripts = list(df_bambu[df_bambu['mgi_symbol'] == gene]['ensembl_transcript_id'])
transcripts.sort()
transcript_color = ['b', 'g']

cond_tx_dfs = []
for this_cond in conditions:
    for this_transcript in transcripts:
        this_filename = os.path.join(data_dir, ds, 'isoforms', gene, f'{this_transcript}_{this_cond}.mAFiA.bed')
        this_df = pd.read_csv(this_filename, sep='\t', dtype={'chrom': str})
        this_df = this_df[this_df['coverage'] >= thresh_coverage]
        this_df.rename(columns={
            'coverage': f'coverage_{this_cond}_{this_transcript}',
            'modRatio': f'modRatio_{this_cond}_{this_transcript}',
            'confidence': f'confidence_{this_cond}_{this_transcript}',
        }, inplace=True)
        cond_tx_dfs.append(this_df)
df_merged = reduce(lambda left, right: pd.merge(left, right, on=bed_fields), cond_tx_dfs)

# df_expanded = []
# for this_cond in conditions:
#     for this_transcript in transcripts:
#         site_mod_ratios = df_merged[f'modRatio_{this_cond}_{this_transcript}'].values
#         df_expanded.append(
#             pd.DataFrame({'modRatio': site_mod_ratios, 'condition': this_cond, 'transcript': this_transcript})
#         )
# df_expanded = pd.concat(df_expanded)
#
# # xticks = [f'{this_cond}\n{this_transcript}' for this_cond in conditions for this_transcript in transcripts]
# yticks = np.arange(0, df_expanded['modRatio'].max(), 25)
#
# plt.figure(figsize=(9*cm, 4*cm))
# sns.stripplot(data=df_expanded, x='condition', y='modRatio', hue='transcript',
#               dodge=True, palette=transcript_color, size=3, rasterized=True)
# plt.legend(title=gene, loc='lower right', bbox_to_anchor=(1.0, 1.0))
# plt.yticks(yticks)
# # plt.title(gene)
# # plt.savefig(os.path.join(img_out, f'strip_plot_isoforms_mod_ratio_{gene}_cov{thresh_coverage}.{FMT}'), **fig_kwargs)

### heatmap ###
df_thresh = df_merged[df_merged[f'modRatio_{conditions[0]}_{transcripts[0]}'] >= 0]
df_thresh.sort_values(f'modRatio_{conditions[0]}_{transcripts[0]}', ascending=False, inplace=True)
mat_isoform_cond = df_thresh[[
    f'modRatio_{this_cond}_{this_transcript}'
for this_transcript in transcripts for this_cond in conditions]].values

# isoform_cond = [', '.join([this_isoform, this_cond]) for this_isoform in transcripts for this_cond in conditions]
cond_isoform = [', '.join([this_cond, this_isoform]) for this_isoform in transcripts for this_cond in conditions]
vec_pos_mod = [', '.join([str(pos), rf'${{{dict_mod_display[mod]}}}$']) for (pos, mod) in df_thresh[['chromEnd', 'name']].values]
chrom = df_thresh['chrom'].unique()[0]

plt.figure(figsize=(6*cm, 4*cm))
ax = plt.gca()
im = ax.imshow(mat_isoform_cond.T, cmap='plasma', vmin=25, vmax=75)
plt.xticks(range(len(vec_pos_mod)), vec_pos_mod, rotation=-90)
plt.yticks(range(4), cond_isoform)
plt.title(f'chr{chrom}, {gene}')
# plt.gca().invert_yaxis()
# plt.gca().xaxis.tick_top()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_ticks([25, 50, 75])
plt.savefig(os.path.join(img_out, f'heatmap_isoforms_mod_ratio_{gene}_cov{thresh_coverage}.{FMT}'), **fig_kwargs)

# for this_row in mat_cond_isoform:
#     plt.plot([0, 1], [this_row[0], this_row[1]], c='b', ls='-')
#     plt.plot([1, 2], [this_row[1], this_row[2]], c='gray', ls='--')
#     plt.plot([2, 3], [this_row[2], this_row[3]], c='b', ls='-')