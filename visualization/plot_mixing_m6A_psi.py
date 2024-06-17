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

mod = 'm6A'
thresh_confidence = 80
thresh_coverage = 50

def import_mAFiA(ds, thresh_conf=80.0, thresh_cov=50):
    df_mAFiA = pd.read_csv(os.path.join(source_data_dir, ds, 'chrALL.mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t')
    df_mAFiA_thresh = df_mAFiA[
        (df_mAFiA['confidence']>=thresh_conf)
        * (df_mAFiA['coverage']>=thresh_cov)
        * (df_mAFiA['name']==mod)
    ]
    return df_mAFiA_thresh


def import_ref(ref_path):
    df_ref = pd.read_csv(ref_path, dtype={'chrom': str}, sep='\t')
    df_ref.rename(columns={'score': 'modRatio'}, inplace=True)
    return df_ref

########################################################################################################################


source_data_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293'
img_out = '/home/adrian/img_out/mixing'
os.makedirs(img_out, exist_ok=True)

all_ds = [
    '0WT',
    '25WT',
    '50WT',
    '75WT',
    # '100WT'
]
ds_names = {this_ds: f"{this_ds.rstrip('WT')}% WT" for this_ds in all_ds}
f_wt = {this_ds: int(this_ds.rstrip('WT'))/100.0 for this_ds in all_ds}

# ref_name = '100% WT'
# df_reference = import_mAFiA('100_WT_0_IVT', thresh_conf=thresh_confidence, thresh_cov=thresh_coverage)

if mod == 'm6A':
    reference_path = '/home/adrian/Data/GLORI/bed_files/GLORI.chrALL.tsv'
    ref_name = 'GLORI'
    markersize = 1
    alpha = 0.5
elif mod == 'psi':
    reference_path = '/home/adrian/Data/PRAISE/PRAISE_HEK293T_span1-3.bed'
    ref_name = 'PRAISE'
    markersize = 3
    alpha = 1.0
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


ticks = np.linspace(0, 100, 5)
xylim = [-1, 101]
fig_scatter, axs = plt.subplots(nrows=1, ncols=len(all_ds), figsize=(6*len(all_ds)*cm, 5*cm))
for ind, this_ds in enumerate(all_ds):
    df_mafia = import_mAFiA(this_ds, thresh_conf=thresh_confidence, thresh_cov=thresh_coverage)
    scatter_x, scatter_y = scatter_mafia_vs_ref(df_reference, df_mafia)
    axs[ind].plot(scatter_x, scatter_y, '.', markersize=markersize, alpha=0.5)
    axs[ind].set_xlim(xylim)
    axs[ind].set_ylim(xylim)
    guideline_x = ticks
    guideline_y = f_wt[this_ds] * guideline_x
    axs[ind].plot(guideline_x, guideline_y, 'r--', alpha=0.5)
    axs[ind].set_xticks(ticks)
    axs[ind].set_yticks(ticks)
    axs[ind].set_xlabel(ref_name)
    # axs[ind].set_ylabel(ds_names[this_ds])
    if ind==0:
        axs[ind].set_ylabel('$\psi$-co-mAFiA')
    axs[ind].set_title(ds_names[this_ds])
# fig_scatter.suptitle(mod, fontsize=15)
fig_scatter.savefig(os.path.join(img_out, f'scatter_mixing_{mod}_conf{thresh_confidence}_cov{thresh_coverage}.{FMT}'), **fig_kwargs)