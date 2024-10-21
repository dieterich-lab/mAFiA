import os
import pandas as pd
import re
import numpy as np
from scipy.stats import kstest
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


def get_logit(vec_modRatio):
    return np.log2(vec_modRatio / (100 - vec_modRatio))


########################################################################################################################
gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
bed_fields = [
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
df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

########################################################################################################################

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

ds_colors = {
    'HFpEF': 'r',
    'Diet': 'b',
    'TAC': 'g'
}

sel_gene = 'Tnnt2'
sel_mod = 'm6A'
sel_gene_row = df_gene[df_gene['gene'] == sel_gene].iloc[0]

# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']

# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

ds = 'Diet'
conditions = ['WT_CD_merged', 'WT_WD_merged']


cond_df = {}
for this_cond in conditions:
    res_file = os.path.join('/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart',
                            f'{ds}/{this_cond}/chrALL.mAFiA.sites.bed')
    df_res = pd.read_csv(res_file, sep='\t', dtype={'chrom': str})
    cond_df[this_cond] = df_res[
        (df_res.chrom == sel_gene_row.chrom)
        * (df_res.chromStart >= sel_gene_row.chromStart)
        * (df_res.chromEnd <= sel_gene_row.chromEnd)
        * (df_res.strand == sel_gene_row.strand)
        * (df_res.name == sel_mod)
        # * (df_res.confidence >= 50.0)
    ]

df_merged = pd.merge(*cond_df.values(),
                     on=['chrom',
                         'chromStart',
                         'chromEnd',
                         'name',
                         'score',
                         'strand',
                         'ref5mer'],
                     suffixes=[f"_{this_cond.rstrip('_merged')}" for this_cond in conditions])

logit1 = get_logit(df_merged[f"modRatio_{conditions[1].rstrip('_merged')}"].values)
logit0 = get_logit(df_merged[f"modRatio_{conditions[0].rstrip('_merged')}"].values)

df_merged['delta_logit'] = logit1 - logit0

vec_delta_logit = df_merged['delta_logit'].values
vec_delta_logit = vec_delta_logit[~np.isinf(vec_delta_logit)]

xlim = [-2.5, 2.5]
bin_width = 0.25
num_bins = int((xlim[1] - xlim[0]) / bin_width)

sigma = np.std(vec_delta_logit)
hist, bin_edges = np.histogram(vec_delta_logit, range=xlim, bins=num_bins)
# norm_hist = hist / hist.sum()
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
vec_x = np.linspace(*xlim, 100)
gaussian = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-vec_x**2 / (2 * sigma**2))
gaussian = gaussian / gaussian.max() * hist.max()

ks_stat, pval = kstest(gaussian, vec_delta_logit)
mean_delta = np.mean(vec_delta_logit)

yticks = np.int64(np.arange(0, np.max(hist), 5))

plt.figure(figsize=(5*cm, 2*cm))
plt.bar(bin_centers, hist, fc=ds_colors[ds])
plt.axvline(x=0, c='gray', ls='--')
plt.plot(vec_x, gaussian, c='gray', ls='--')
plt.xlim(xlim)
plt.yticks(yticks)
plt.text(0.01, 0.95, f'mean: {mean_delta:.3f}\np-val: {pval:.3E}',
         ha='left', va='top', transform=plt.gca().transAxes)
plt.savefig(os.path.join(img_out, f"delta_logit_{sel_gene}_{sel_mod}_{ds}.{FMT}"), **fig_kwargs)

