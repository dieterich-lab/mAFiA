import os
HOME = os.path.expanduser('~')
from glob import glob
import pandas as pd
import numpy as np
from scipy.stats import linregress
#######################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt

cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

def filter_df(in_df, sel_motif):
    df_filtered = in_df[
        (np.float32(in_df['P_adjust']) <= P_VAL_THRESH)
        & (in_df['ref_motif']==sel_motif)
        # & (df['pred_motif']==sel_motif)
    ]
    return df_filtered

def get_norm_counts(in_df):
    num_bins = 100
    bin_edges = np.linspace(0, 1, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_counts, _ = np.histogram(in_df['mod_prob'], bins=bin_edges)
    norm_counts = bin_counts / np.sum(bin_counts)

    return norm_counts, bin_centers

def calc_mod_ratio_with_margin(in_df_motif, prob_margin, thresh_cov=50):
    df_motif_avg = pd.DataFrame()
    for index in in_df_motif['index'].unique():
        df_site = in_df_motif[in_df_motif['index']==index]
        df_site_margin = df_site[
            (df_site['mod_prob']<prob_margin)
            | (df_site['mod_prob']>=(1-prob_margin))
            ]
        if len(df_site_margin)<thresh_cov:
            continue
        num_features = len(df_site_margin)
        mod_ratio = np.mean((df_site_margin['mod_prob'].values)>=(1-prob_margin))
        new_row = df_site_margin.iloc[0][:16]
        new_row['num_features'] = num_features
        new_row['mod_ratio'] = mod_ratio
        df_motif_avg = pd.concat([df_motif_avg, new_row.to_frame().T])
    df_motif_avg.reset_index(drop=True)
    return df_motif_avg

# train_dataset = 'ISA-WUE'
# results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'

train_dataset = 'm6Anet'
results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6Anet'

test_datasets = [
    '0_WT_100_IVT',
    '25_WT_75_IVT',
    '50_WT_50_IVT',
    '75_WT_25_IVT',
    '100_WT_0_IVT'
    # 'P2_WT'
]
ds_colors = {
    '0_WT_100_IVT' : 'b',
    '25_WT_75_IVT' : 'g',
    '50_WT_50_IVT' : 'm',
    '75_WT_25_IVT' : 'c',
    '100_WT_0_IVT' : 'r',
    'P2_WT': 'r'
}

ds_names = {
    '0_WT_100_IVT' : '0% WT',
    '25_WT_75_IVT' : '25% WT',
    '50_WT_50_IVT' : '50% WT',
    '75_WT_25_IVT' : '75% WT',
    '100_WT_0_IVT' : '100% WT',
    'P2_WT' : '100% WT P2'
}

img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('mRNA_train_{}_test_HEK293_GLORI'.format(train_dataset)))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

P_VAL_THRESH = 1.0E-99
COV_THRESH = 50
PROB_MARGIN = 0.5
COMMON_SITES_ONLY = False

def import_MAFIA_res(res_dir):
    dfs = {}
    for ds in test_datasets:
        paths = glob(os.path.join(res_dir, 'res_train_{}_test_{}.tsv.merged'.format(train_dataset, ds)))
        if len(paths)==1:
            df = pd.read_csv(paths[0], sep='\t').rename(columns={'Unnamed: 0': 'index'})
            dfs[ds] = df
        elif len(paths)>=2:
            df = pd.concat([pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'}) for path in paths])
            dfs[ds] = df
    return dfs

def import_m6Anet_res(res_dir):
    dfs = {}
    for ds in test_datasets:
        path = os.path.join(res_dir, ds, 'data.site_proba.glori_filtered.csv')
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        df = df.rename(columns={"kmer": "ref_motif", "Pvalue": "P_adjust"})
        dfs[ds] = df
    return dfs


# motifs = list(set.intersection(*[set(df['ref_motif'].unique()) for df in dfs.values()]))
motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']

### aggregate mod. ratio per site ###
# def calc_mod_ratio(in_df_motif, thresh_mod=0.5, thresh_cov=50):
#     df_motif_avg = pd.DataFrame()
#     for index in in_df_motif['index'].unique():
#         df_site = in_df_motif[in_df_motif['index']==index]
#         if len(df_site)<thresh_cov:
#             continue
#         num_features = len(df_site)
#         mod_ratio = np.mean((df_site['mod_prob'].values)>=thresh_mod)
#         new_row = df_site.iloc[0][:16]
#         new_row['num_features'] = num_features
#         new_row['mod_ratio'] = mod_ratio
#         df_motif_avg = pd.concat([df_motif_avg, new_row.to_frame().T])
#     df_motif_avg.reset_index(drop=True)
#     return df_motif_avg
#
# df_motif_avg_ivt = calc_mod_ratio(df_motif_ivt, thresh_mod=crit_thresh, thresh_cov=COV_THRESH)
# df_motif_avg_wt = calc_mod_ratio(df_motif_wt, thresh_mod=crit_thresh, thresh_cov=COV_THRESH)

# dfs = import_MAFIA_res(results_dir)
dfs = import_m6Anet_res(results_dir)

### plots ###
num_rows = 3
num_cols = 2
# fig_width = num_cols*5
# fig_height = num_rows*5
fig_width = 6*cm
fig_height = 10*cm
xticks = [0, 0.5, 1]
yticks = xticks
# fig_hist = plt.figure(figsize=(fig_width, fig_height))
# axes_hist = fig_hist.subplots(num_rows, num_cols).flatten()
fig_mod_ratio = plt.figure(figsize=(fig_width, fig_height))
axes_mod_ratio = fig_mod_ratio.subplots(num_rows, num_cols).flatten()
dict_ds_motif_avg = {}
for ds, df in dfs.items():
    dict_ds_motif_avg[ds] = {}
    for subplot_ind, this_motif in enumerate(motifs):
        if this_motif not in df['ref_motif'].unique():
            continue
        ds_motif = filter_df(df, this_motif)
        # axes_hist[subplot_ind].step(ds_bin_centers, ds_norm_counts, color=ds_colors[ds], label='{}'.format(ds))
        #
        # axes_hist[subplot_ind].axvspan(xmin=PROB_MARGIN, xmax=1 - PROB_MARGIN, color='gray', alpha=0.5)
        # if subplot_ind>=num_cols*(num_rows-1):
        #     axes_hist[subplot_ind].set_xlabel('Mod. Prob.')
        # if subplot_ind%num_cols==0:
        #     axes_hist[subplot_ind].set_ylabel('Norm. frequency')
        # axes_hist[subplot_ind].set_title('{}'.format(this_motif))
        # axes_hist[subplot_ind].set_xlim([-0.05, 1.05])
        # axes_hist[subplot_ind].set_ylim([0, 0.3])
        # axes_hist[subplot_ind].legend(loc='upper left')

        if train_dataset not in ['m6Anet', 'CHEUI']:
            ds_norm_counts, ds_bin_centers = get_norm_counts(df)
            df_motif_avg = calc_mod_ratio_with_margin(ds_motif, prob_margin=PROB_MARGIN, thresh_cov=COV_THRESH)
        else:
            df_motif_avg = ds_motif
        dict_ds_motif_avg[ds][this_motif] = df_motif_avg

        glori_ratio = np.float64(df_motif_avg['Ratio'])
        mod_ratio = np.float64(df_motif_avg['mod_ratio'])
        corr = np.corrcoef(glori_ratio, mod_ratio)[0, 1]

        if ds in ['100_WT_0_IVT', 'P2_WT']:
            axes_mod_ratio[subplot_ind].scatter(glori_ratio, mod_ratio, color=ds_colors[ds], marker='.', s=1, label=f'corr. {corr:.2f}')
            # axes_mod_ratio[subplot_ind].scatter(glori_ratio, mod_ratio, color=ds_colors[ds], marker='.', label='{} sites, corr. {:.2f}'.format(len(glori_ratio), corr))
            # axes_mod_ratio[subplot_ind].plot(glori_ratio_ivt, mod_ratio_ivt, 'b.', label='IVT, {} sites, corr. {:.2f}'.format(len(glori_ratio_ivt), corr_ivt))
            # axes_mod_ratio[subplot_ind].plot(glori_ratio_wt, mod_ratio_wt, 'r.', label='WT, {} sites, corr. {:.2f}'.format(len(glori_ratio_wt), corr_wt))

            axes_mod_ratio[subplot_ind].plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'k--', alpha=0.5)
            axes_mod_ratio[subplot_ind].set_xlim([-0.05, 1.05])
            axes_mod_ratio[subplot_ind].set_ylim([-0.05, 1.05])
            if subplot_ind >= num_cols * (num_rows - 1):
                # axes_mod_ratio[subplot_ind].set_xlabel('GLORI mod. ratio')
                axes_mod_ratio[subplot_ind].set_xticks(xticks)
            else:
                axes_mod_ratio[subplot_ind].set_xticks([])
            if subplot_ind%num_cols==0:
                # axes_mod_ratio[subplot_ind].set_ylabel('ONT mod. ratio')
                axes_mod_ratio[subplot_ind].set_yticks(yticks)
            else:
                axes_mod_ratio[subplot_ind].set_yticks([])
            axes_mod_ratio[subplot_ind].set_title(f'{this_motif}', pad=-10)
            axes_mod_ratio[subplot_ind].legend(loc='upper left')

            print(this_motif, len(glori_ratio))

# fig_hist.tight_layout()
# fig_hist.savefig(os.path.join(img_out, f'hist_modProbs_pValThresh{P_VAL_THRESH}_marginProb{PROB_MARGIN}.{FMT}'), **fig_kwargs)

fig_mod_ratio.tight_layout()
# fig_mod_ratio.suptitle('Train: {}\nTest: {}'.format(train_dataset, 'HEK293 WT'))
fig_mod_ratio.savefig(os.path.join(img_out, f'corr_glori_modRatio_pValThresh{P_VAL_THRESH}_covThresh{COV_THRESH}_marginProb{PROB_MARGIN}.{FMT}'), **fig_kwargs)

### 100 WT overall ###
df_wt = pd.concat(dict_ds_motif_avg['100_WT_0_IVT'].values())
corr = np.corrcoef(df_wt['Ratio'], df_wt['mod_ratio'])[0, 1]
num_sites = len(df_wt)
with open(os.path.join(img_out, 'stats.txt'), 'w') as out_f:
    out_f.write(f'{num_sites} sites, corr {corr:.3f}')

### 2D density plots ###
ds_slope_err = {}
num_bins = 20
interval = 100 / num_bins
fig = plt.figure(figsize=(20*cm, 5*cm))
for subplot_ind, ds in enumerate(dict_ds_motif_avg.keys()):
    ax = plt.subplot(1, len(test_datasets), subplot_ind+1)
    ds_agg_avg = pd.concat(dict_ds_motif_avg[ds].values())
    glori_ratio = np.float32(ds_agg_avg['Ratio'].values)
    ont_ratio = np.float32(ds_agg_avg['mod_ratio'].values)
    hist, x_bins, y_bins = np.histogram2d(glori_ratio, ont_ratio, bins=num_bins, range=[[0, 1], [0, 1]], density=True)
    im = ax.imshow(hist.T, origin='lower', vmin=1, vmax=6, cmap='plasma')
    ax.set_xticks(np.arange(0, num_bins+1, 5)-0.5, np.int32(np.arange(0, num_bins+1, 5)*interval))
    for tick in ax.yaxis.get_majorticklabels():
        tick.set_verticalalignment('bottom')
    # ax.set_xlabel('GLORI')
    if subplot_ind==0:
        # ax.set_ylabel('ONT mod. ratio')
        ax.set_yticks(np.arange(0, num_bins + 1, 5) - 0.5, np.int32(np.arange(0, num_bins + 1, 5) * interval))
    else:
        ax.set_yticks([])

    # ax.set_title('{}\n{} sites'.format(ds, len(ds_agg_avg)))
    ax.set_title(ds_names[ds])

    lin_fit = linregress(0.5 * (x_bins[1:] + x_bins[:-1]), x_bins[np.argmax(hist, axis=1)])
    ds_slope_err[ds] = (lin_fit.slope, lin_fit.stderr)
# fig.tight_layout()
cb_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
cbar = fig.colorbar(im, cax=cb_ax)
# cb_ax.set_ylabel('Norm. site count', rotation=-90, labelpad=20)
fig.savefig(os.path.join(img_out, f'hist2d_glori_modRatio_pValThresh{P_VAL_THRESH}_covThresh{COV_THRESH}_marginProb{PROB_MARGIN}.{FMT}'), **fig_kwargs)

### slope vs mixing ratio ###
ticks = np.arange(0, 1.01, 0.25)
plt.figure(figsize=(5*cm, 5*cm))
plt.errorbar(x=ticks,
             y=[v[0] for v in ds_slope_err.values()],
             yerr=[v[1] for v in ds_slope_err.values()],
             marker='.',
             linestyle='None',
             capsize=2.0
             )
plt.xticks(ticks)
plt.yticks(ticks)
# plt.xlabel('WT Ratio')
# plt.ylabel('Ridge Inclination')
plt.savefig(os.path.join(img_out, f'ridge_inclination_pValThresh{P_VAL_THRESH}_covThresh{COV_THRESH}_marginProb{PROB_MARGIN}.{FMT}'), **fig_kwargs)

# plt.close('all')