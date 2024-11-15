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
import seaborn as sns
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
sns.set_style('whitegrid')

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

data_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/mouse_heart'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'

ds = 'HFpEF'
conditions = ['ctrl', 'HFpEF']
bambu_file = os.path.join(
    '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/DTU',
    'DTU_HFpEF_Laura_ENS102.tsv'
)

df_bambu = pd.read_csv(bambu_file, sep='\t')

gene = 'Mfn2'

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

df_expanded = []
for this_cond in conditions:
    for this_transcript in transcripts:
        site_mod_ratios = df_merged[f'modRatio_{this_cond}_{this_transcript}'].values
        df_expanded.append(
            pd.DataFrame({'modRatio': site_mod_ratios, 'condition': this_cond, 'transcript': this_transcript})
        )
df_expanded = pd.concat(df_expanded)

# xticks = [f'{this_cond}\n{this_transcript}' for this_cond in conditions for this_transcript in transcripts]
yticks = np.arange(0, df_expanded['modRatio'].max(), 25)

plt.figure(figsize=(9*cm, 4*cm))
sns.stripplot(data=df_expanded, x='condition', y='modRatio', hue='transcript',
              dodge=True, palette=transcript_color, size=3, rasterized=True)
plt.legend(title=gene, loc='lower right', bbox_to_anchor=(1.0, 1.0))
plt.yticks(yticks)
# plt.title(gene)
# plt.savefig(os.path.join(img_out, f'strip_plot_isoforms_mod_ratio_{gene}_cov{thresh_coverage}.{FMT}'), **fig_kwargs)