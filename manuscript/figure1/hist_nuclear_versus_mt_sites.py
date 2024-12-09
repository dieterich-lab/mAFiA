import os
import pandas as pd
import numpy as np
from tqdm import tqdm
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

sns.set_style('whitegrid')

mod_color = {
    'm6A': 'r',
    'psi': 'b'
}
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

gene_bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'frame',
    'description'
]

protein_gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/protein_coding.gene.GRCm38.102.bed'
df_protein_gene = pd.read_csv(protein_gene_bed, sep='\t', names=gene_bed_fields)

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

thresh_conf = 80
thresh_cov = 50

ds = 'TAC'
condition = 'SHAM_merged'
# ds = 'HFpEF'
# condition = 'HFpEF_merged'

in_filename = os.path.join(res_dir, ds, condition, 'chrALL.mAFiA.sites.bed')
out_filename = os.path.join(res_dir, ds, condition, f'protein_coding.conf{thresh_conf}_cov{thresh_cov}.chrALL.mAFiA.sites.bed')

if os.path.exists(out_filename):
    df_site_thresh_annotated = pd.read_csv(out_filename, sep='\t', dtype={'chrom': str})
else:
    df_sites = pd.read_csv(in_filename, sep='\t', dtype={'chrom': str})
    df_site_thresh = df_sites[
        (df_sites['coverage'] >= thresh_cov)
        * (df_sites['confidence'] >= thresh_conf)
    ]

    df_site_thresh['annotation'] = None
    for _, this_row in tqdm(df_protein_gene.iterrows()):
        mask = (df_site_thresh['chrom'] == this_row['chrom'])\
               * (df_site_thresh['chromStart'] >= this_row['chromStart'])\
               * (df_site_thresh['chromEnd'] <= this_row['chromEnd'])\
               * (df_site_thresh['strand'] == this_row['strand'])
        df_site_thresh.loc[mask, 'annotation'] = this_row['name']

    df_site_thresh_annotated = df_site_thresh[~df_site_thresh['annotation'].isna()]
    df_site_thresh_annotated.to_csv(out_filename, sep='\t', index=False)

# df_site_thresh_nuc = df_site_thresh[
#     df_site_thresh['chrom'].isin([str(this_chr) for this_chr in (range(1, 20))] + ['X', 'Y'])
# ]
# hist_nuc, _ = np.histogram(df_site_thresh_nuc['modRatio'].values, range=[0, 100], bins=100)
# df_site_thresh_mt = df_site_thresh[
#     df_site_thresh['chrom'] == 'MT'
# ]

# plt.figure(figsize=(5*cm, 5*cm))
# plt.hist(df_site_thresh_nuc['modRatio'], , fc='b', alpha=0.5)
# plt.hist(df_site_thresh_mt['modRatio'], range=[0, 100], bins=100, fc='r', alpha=0.5)

plt.figure(figsize=(12*cm, 4*cm))
sns.stripplot(x='chrom', y='modRatio', data=df_site_thresh_annotated,
              hue='name', hue_order=list(dict_mod_display.keys()),
              palette=mod_color, dodge=True, jitter=0.2, size=1, rasterized=True)
mod_labels = [rf'${{{this_val}}}$' for this_val in dict_mod_display.values()]
plt.ylim([0, 100])
plt.yticks(np.linspace(0, 100, 5))
plt.legend(title=None, labels=mod_labels, loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(img_out, f'strip_plot_sites_by_chrom_{ds}_{condition}.{FMT}'), **fig_kwargs)
