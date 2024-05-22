import os
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.stats import linregress

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'pdf'
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################


def import_mAFiA(ds, thresh_conf=80.0, thresh_cov=50):
    df_mAFiA = pd.read_csv(os.path.join(source_data_dir, ds, 'chrALL.mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t')
    df_mAFiA_thresh = df_mAFiA[
        (df_mAFiA['confidence']>=thresh_conf)
        * (df_mAFiA['coverage']>=thresh_cov)
    ]
    return df_mAFiA_thresh


def import_ref(ref_path):
    df_ref = pd.read_csv(ref_path, dtype={'chrom': str}, sep='\t')
    df_ref.rename(columns={'score': 'modRatio'}, inplace=True)
    # df_glori['chrom'] = [chr.lstrip('chr') for chr in df_glori['Chr']]
    # df_glori['chromStart'] = df_glori['Sites'] - 1
    # df_glori['chromEnd'] = df_glori['Sites']
    # df_glori['strand'] = df_glori['Strand']
    # df_glori['modRatio'] = np.int32(np.round(df_glori['NormeRatio'] * 100.0))
    # df_glori_thresh = df_glori[df_glori['P_adjust']<thresh_pval]

    return df_ref

########################################################################################################################


source_data_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293'
img_out = '/home/adrian/img_out/psi_mixing'
os.makedirs(img_out, exist_ok=True)

all_ds = [
    '0_WT_100_IVT',
    '25_WT_75_IVT',
    '50_WT_50_IVT',
    '75_WT_25_IVT',
    '100_WT_0_IVT'
]
ds_names = {this_ds: f"{this_ds.split('_')[0]}% WT" for this_ds in all_ds}
f_wt = {this_ds: int(this_ds.split('_')[0])/100.0 for this_ds in all_ds}

# reference_path = '/home/adrian/Data/PRAISE/PRAISE_HEK293T.bed'
reference_path = '/home/adrian/Data/PRAISE/PRAISE_HEK293T_span1-3.bed'
ref_name = 'PRAISE'
# reference_path = '/home/adrian/Data/BID_seq/BID_seq_HEK293T.bed'
# ref_name = 'BID-Seq'
df_reference = import_ref(reference_path)

########################################################################################################################
### Density plots ######################################################################################################
########################################################################################################################
def hist2d_mafia_vs_ref(df_ref, df_pred, n_bins, drop_margin):
    df_merged = pd.merge(df_ref, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=[f'_{ref_name}', '_mafia'])
    counts, bin_y, bin_x = np.histogram2d(
        df_merged['modRatio_mafia'], df_merged[f'modRatio_{ref_name}'],
        bins=[n_bins, n_bins], range=[[0, 100], [0, 100]],
    )
    log_counts = np.log10(1 + counts)

    fit_vec_x = 0.5 * (bin_x[1:] + bin_x[:-1])
    fit_vec_y = fit_vec_x[np.argmax(counts, axis=0)]
    lin_fit = linregress(fit_vec_x[drop_margin:-drop_margin], fit_vec_y[drop_margin:-drop_margin])

    return counts, log_counts, bin_x, bin_y, lin_fit.slope, lin_fit.stderr

def scatter_mafia_vs_ref(df_ref, df_pred):
    df_merged = pd.merge(df_ref, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=[f'_{ref_name}', '_mafia'])
    return df_merged[f'modRatio_{ref_name}'].values, df_merged['modRatio_mafia'].values


thresh_confidence = 80
thresh_coverage = 100
ticks = np.linspace(0, 100, 5)
xylim = [-1, 101]
fig_scatter, axs = plt.subplots(nrows=1, ncols=len(all_ds), figsize=(6*len(all_ds)*cm, 5*cm))
for ind, this_ds in enumerate(all_ds):
    df_mafia = import_mAFiA(this_ds, thresh_conf=thresh_confidence, thresh_cov=thresh_coverage)
    scatter_x, scatter_y = scatter_mafia_vs_ref(df_reference, df_mafia)
    axs[ind].plot(scatter_x, scatter_y, '.')
    axs[ind].set_xlim(xylim)
    axs[ind].set_ylim(xylim)
    guideline_x = ticks
    guideline_y = f_wt[this_ds] * guideline_x
    axs[ind].plot(guideline_x, guideline_y, 'r--', alpha=0.5)
    axs[ind].set_xticks(ticks)
    axs[ind].set_yticks(ticks)
    axs[ind].set_xlabel(ref_name)
    if ind==0:
        axs[ind].set_ylabel('$\psi$-co-mAFiA')
    axs[ind].set_title(ds_names[this_ds])
fig_scatter.savefig(os.path.join(img_out, f'scatter_mixing_mAFiA_vs_{ref_name}_conf{thresh_confidence}_cov{thresh_coverage}.{FMT}'), **fig_kwargs)