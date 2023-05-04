import os
HOME = os.path.expanduser('~')
from glob import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def get_norm_counts(in_df, sel_motif):
    df_motif = in_df[
        (in_df['P_adjust'] <= P_VAL_THRESH)
        & (in_df['ref_motif']==sel_motif)
        # & (df['pred_motif']==sel_motif)
    ]

    ### histogram of mod prob per read ###
    num_bins = 100
    bin_edges = np.linspace(0, 1, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_counts, _ = np.histogram(df_motif['mod_prob'], bins=bin_edges)
    norm_counts = bin_counts / np.sum(bin_counts)

    return norm_counts, bin_centers, df_motif

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


# parser = argparse.ArgumentParser()
# parser.add_argument('--df_file')
# args = parser.parse_args()
# df_file = args.df_file

# classifier = 'random_ligation_A_m6A_MaxAbs'
# ivt_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/old/results/res_HEK293_IVT_{}_noEnforceMotif_modProbPerRead.tsv'.format(classifier)
# wt_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_{}_modProbThresh0.5.tsv'.format(classifier)
# ivt_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_0_WT_100_IVT_RTA_20230221_WUE_splint_lig_modProbPerRead.tsv.merged'
# wt_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_100_WT_0_IVT_RTA_20230221_WUE_splint_lig_modProbPerRead.tsv.merged'
# img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('HEK293_WT_vs_IVT_{}'.format(classifier)))

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
training_dataset = 'WUE_splint_batch*'
testing_datasets = [
    '0_WT_100_IVT_RTA',
    '25_WT_75_IVT_RTA',
    '50_WT_50_IVT_RTA',
    '75_WT_25_IVT_RTA',
    '100_WT_0_IVT_RTA'
]
ds_colors = {
    '0_WT_100_IVT_RTA' : 'b',
    '25_WT_75_IVT_RTA' : 'g',
    '50_WT_50_IVT_RTA' : 'm',
    '75_WT_25_IVT_RTA' : 'c',
    '100_WT_0_IVT_RTA' : 'r'
}

img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('HEK293_mixing_{}'.format(training_dataset)))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

P_VAL_THRESH = 1.0E-99
COV_THRESH = 50
PROB_MARGIN = 0.2
COMMON_SITES_ONLY = False

# dfs = {
#     ds : pd.read_csv(os.path.join(results_dir, 'res_{}_{}_modProbPerRead.tsv.merged'.format(ds, training_dataset)), sep='\t').rename(columns={'Unnamed: 0': 'index'}) for ds in testing_datasets
# }
dfs = {}
for ds in testing_datasets:
    paths = glob(os.path.join(results_dir, 'res_{}_{}_modProbPerRead.tsv.merged'.format(ds, training_dataset)))
    if len(paths)==1:
        df = pd.read_csv(paths[0], sep='\t').rename(columns={'Unnamed: 0': 'index'})
        dfs[ds] = df
    elif len(paths)>=2:
        df = pd.concat([pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'}) for path in paths])
        dfs[ds] = df

# motifs = list(set.intersection(*[set(df['ref_motif'].unique()) for df in dfs.values()]))
motifs = ['AGACT', 'GGACA', 'GGACC', 'GAACT', 'GGACT', 'TGACT']

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


### plots ###
num_rows = 2
num_cols = 3
fig_width = 16
fig_height = 10
fig_hist = plt.figure(figsize=(fig_width, fig_height))
axes_hist = fig_hist.subplots(num_rows, num_cols).flatten()
fig_mod_ratio = plt.figure(figsize=(fig_width, fig_height))
axes_mod_ratio = fig_mod_ratio.subplots(num_rows, num_cols).flatten()
dict_ds_motif_avg = {}
for ds, df in dfs.items():
    dict_ds_motif_avg[ds] = {}
    for subplot_ind, this_motif in enumerate(motifs):
        if this_motif not in df['ref_motif'].unique():
            continue
        ds_norm_counts, ds_bin_centers, ds_motif = get_norm_counts(df, this_motif)
        axes_hist[subplot_ind].step(ds_bin_centers, ds_norm_counts, color=ds_colors[ds], label='{}, {} features'.format(ds, len(ds_motif)))

        axes_hist[subplot_ind].axvspan(xmin=PROB_MARGIN, xmax=1 - PROB_MARGIN, color='gray', alpha=0.5)
        axes_hist[subplot_ind].set_xlabel('Mod. Prob.', fontsize=15)
        axes_hist[subplot_ind].set_ylabel('Norm. frequency', fontsize=15)
        axes_hist[subplot_ind].set_title('{}'.format(this_motif), fontsize=20)
        axes_hist[subplot_ind].set_xlim([-0.05, 1.05])
        axes_hist[subplot_ind].set_ylim([0, 0.3])
        axes_hist[subplot_ind].legend(loc='upper left', fontsize=12)

        df_motif_avg = calc_mod_ratio_with_margin(ds_motif, prob_margin=PROB_MARGIN, thresh_cov=COV_THRESH)
        dict_ds_motif_avg[ds][this_motif] = df_motif_avg

        glori_ratio = np.float64(df_motif_avg['Ratio'])
        mod_ratio = np.float64(df_motif_avg['mod_ratio'])
        corr = np.corrcoef(glori_ratio, mod_ratio)[0, 1]

        if ds=='100_WT_0_IVT_RTA':
            axes_mod_ratio[subplot_ind].scatter(glori_ratio, mod_ratio, color=ds_colors[ds], marker='.', label='{}, {} sites, corr. {:.2f}'.format(ds, len(glori_ratio), corr))
            # axes_mod_ratio[subplot_ind].plot(glori_ratio_ivt, mod_ratio_ivt, 'b.', label='IVT, {} sites, corr. {:.2f}'.format(len(glori_ratio_ivt), corr_ivt))
            # axes_mod_ratio[subplot_ind].plot(glori_ratio_wt, mod_ratio_wt, 'r.', label='WT, {} sites, corr. {:.2f}'.format(len(glori_ratio_wt), corr_wt))

            axes_mod_ratio[subplot_ind].plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'k--', alpha=0.5)
            axes_mod_ratio[subplot_ind].set_xlim([-0.05, 1.05])
            axes_mod_ratio[subplot_ind].set_ylim([-0.05, 1.05])
            axes_mod_ratio[subplot_ind].set_xlabel('GLORI mod. ratio', fontsize=15)
            axes_mod_ratio[subplot_ind].set_ylabel('ONT mod. ratio', fontsize=15)
            axes_mod_ratio[subplot_ind].set_title('{}'.format(this_motif), fontsize=20)
            axes_mod_ratio[subplot_ind].legend(loc='upper left', fontsize=12)
fig_hist.tight_layout()
fig_mod_ratio.tight_layout()

fig_hist.savefig(os.path.join(img_out, 'hist_modProbs_pValThresh{:.2E}_marginProb{:.2f}.png'.format(P_VAL_THRESH, PROB_MARGIN)), bbox_inches='tight')
fig_mod_ratio.savefig(os.path.join(img_out, 'corr_glori_modRatio_pValThresh{:.2E}_covThresh{}_marginProb{:.2f}.png'.format(P_VAL_THRESH, COV_THRESH, PROB_MARGIN)), bbox_inches='tight')

### 2D density plots ###
num_bins = 20
interval = 100 / num_bins
fig = plt.figure(figsize=(len(testing_datasets)*5, 5))
for subplot_ind, ds in enumerate(dict_ds_motif_avg.keys()):
    ax = plt.subplot(1, len(testing_datasets), subplot_ind+1)
    ds_agg_avg = pd.concat(dict_ds_motif_avg[ds].values())
    glori_ratio = np.float32(ds_agg_avg['Ratio'].values)
    ont_ratio = np.float32(ds_agg_avg['mod_ratio'].values)
    hist, x_bins, y_bins = np.histogram2d(glori_ratio, ont_ratio, bins=num_bins, range=[[0, 1], [0, 1]], density=True)
    im = ax.imshow(hist.T, origin='lower', vmin=0.5, vmax=5, cmap='plasma')
    ax.set_xticks(np.arange(0, num_bins+1, 5)-0.5, np.int32(np.arange(0, num_bins+1, 5)*interval))
    ax.set_yticks(np.arange(0, num_bins+1, 5)-0.5, np.int32(np.arange(0, num_bins+1, 5)*interval))
    for tick in ax.yaxis.get_majorticklabels():
        tick.set_verticalalignment('bottom')
    ax.set_xlabel('GLORI', fontsize=15)
    if subplot_ind==0:
        ax.set_ylabel('ONT mod. ratio', fontsize=15)
    ax.set_title('{}\n{} sites'.format(ds, len(ds_agg_avg)), fontsize=15)
# fig.tight_layout()
cb_ax = fig.add_axes([0.92, 0.16, 0.02, 0.68])
cbar = fig.colorbar(im, cax=cb_ax)
cb_ax.set_ylabel('Norm. site count', fontsize=15, rotation=-90, labelpad=20)
fig.savefig(os.path.join(img_out, 'hist2d_glori_modRatio_pValThresh{:.2E}_covThresh{}_marginProb{:.2f}.png'.format(P_VAL_THRESH, COV_THRESH, PROB_MARGIN)), bbox_inches='tight')

# plt.close('all')