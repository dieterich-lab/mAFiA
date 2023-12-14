import os
import pandas as pd
import numpy as np
import matplotlib as mpl
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
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################


def import_mAFiA(ds, thresh_coverage=10):
    df_mAFiA = pd.read_csv(os.path.join(source_data_dir, f'{ds}.mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t')
    df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage']>=thresh_coverage]
    return df_mAFiA_thresh


def import_glori(thresh_pval=1.0):
    df_glori = pd.read_csv(glori_path, dtype={'Chr': str})
    df_glori['chrom'] = [chr.lstrip('chr') for chr in df_glori['Chr']]
    df_glori['chromStart'] = df_glori['Sites'] - 1
    df_glori['chromEnd'] = df_glori['Sites']
    df_glori['strand'] = df_glori['Strand']
    df_glori['modRatio'] = np.int32(np.round(df_glori['NormeRatio'] * 100.0))
    df_glori_thresh = df_glori[df_glori['P_adjust']<thresh_pval]

    return df_glori_thresh

########################################################################################################################


source_data_dir = '/home/adrian/NCOMMS_revision/source_data/MIXING'
img_out = '/home/adrian/NCOMMS_revision/images/MIXING'
os.makedirs(img_out, exist_ok=True)

glori_path = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.ref5mers.csv'

all_ds = [
    '0WT',
    '25WT',
    '50WT',
    '75WT',
    '100WT'
]

f_wt = {this_ds: int(this_ds.rstrip('WT'))/100.0 for this_ds in all_ds}

df_glori = import_glori()

########################################################################################################################
### Density plots ######################################################################################################
########################################################################################################################
def hist2d_mafia_vs_glori(df_ref, df_pred, n_bins, drop_margin):
    df_merged = pd.merge(df_ref, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'strand', 'ref5mer'], suffixes=['_glori', '_mafia'])
    counts, bin_y, bin_x = np.histogram2d(
        df_merged['modRatio_mafia'], df_merged['modRatio_glori'],
        bins=[n_bins, n_bins], range=[[0, 100], [0, 100]],
    )
    log_counts = np.log10(1 + counts)

    fit_vec_x = 0.5 * (bin_x[1:] + bin_x[:-1])
    fit_vec_y = fit_vec_x[np.argmax(counts, axis=0)]
    lin_fit = linregress(fit_vec_x[drop_margin:-drop_margin], fit_vec_y[drop_margin:-drop_margin])

    return counts, log_counts, bin_x, bin_y, lin_fit.slope, lin_fit.stderr

    # df_merged_fit = df_merged[df_merged['modRatio_glori']>=50]
    # df_merged_fit = df_merged
    # def fit_func(x, A):  # this is your 'straight line' y=f(x)
    #     return A * x
    # popt, pcov = curve_fit(fit_func, df_merged_fit['modRatio_glori'], df_merged_fit['modRatio_mafia'])
    # return counts, bin_x, bin_y, popt, pcov

vmax = 50
# vmax = 2
num_bins = 20
margin = 1
ticks = np.linspace(0, num_bins, 5)
ticklabels = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)

# vec_x = np.arange(0, num_bins)

ds_fits = {}
fig_hist2d, axs = plt.subplots(nrows=1, ncols=5, figsize=(20*cm, 4*cm))
for ind, this_ds in enumerate(all_ds):
    df_mafia = import_mAFiA(this_ds)
    this_counts, this_log_counts, this_bin_x, this_bin_y, this_slope, this_stderr = hist2d_mafia_vs_glori(df_glori, df_mafia, num_bins, margin)
    ds_fits[this_ds] = (this_slope, this_stderr)
    im = axs[ind].imshow(this_counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
    # im = axs[ind].imshow(this_log_counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)

    x_vec = np.array([0, np.where(this_bin_x==100)[0][0]])
    y_vec = np.array([0, np.where(this_bin_y==f_wt[this_ds]*100)[0][0]])
    if this_ds=='0WT':
        axs[ind].plot(x_vec-0.5, y_vec, c='r', linestyle='--', alpha=0.5)
    else:
        axs[ind].plot(x_vec-0.5, y_vec-0.5, c='r', linestyle='--', alpha=0.5)

    axs[ind].set_xticks(ticks-0.5, ticklabels)
    axs[ind].set_yticks(ticks-0.5, ticklabels)

    # if ind==0:
    #     axs[ind].set_yticks(ticks-0.5, ticklabels)
    # else:
    #     axs[ind].set_yticks(ticks-0.5, [])

    # plt.plot(vec_x, this_popt[0]*vec_x, c='r', linestyle='--', alpha=0.5)
cbar = fig_hist2d.colorbar(im,  ax=axs, orientation='vertical', location='right', fraction=0.015, aspect=10, pad=0.03)
cbar.set_ticks(np.linspace(0, vmax, 3))
fig_hist2d.savefig(os.path.join(img_out, f'hist2d_mixing_mAFiA_vs_GLORI.{FMT}'), **fig_kwargs)


########################################################################################################################
### Slope versus concentration #########################################################################################
########################################################################################################################
conc_ticks = np.arange(0, 1.01, 0.25)
plt.figure(figsize=(5*cm, 5*cm))
plt.errorbar(x=conc_ticks,
             y=[v[0] for v in ds_fits.values()],
             yerr=[v[1] for v in ds_fits.values()],
             marker='.',
             linestyle='None',
             capsize=2.0
             )
plt.xticks(conc_ticks)
plt.yticks(conc_ticks)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.savefig(os.path.join(img_out, f'slope_vs_concentration.{FMT}'), **fig_kwargs)