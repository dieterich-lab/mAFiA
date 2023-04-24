import os
HOME = os.path.expanduser('~')
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# parser = argparse.ArgumentParser()
# parser.add_argument('--df_file')
# args = parser.parse_args()
# df_file = args.df_file

classifier = 'random_ligation_A_m6A'
ivt_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293_IVT_{}_modProbThresh0.5.tsv'.format(classifier)
wt_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_{}_modProbThresh0.5.tsv'.format(classifier)

img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('HEK293_WT_vs_IVT_{}'.format(classifier)))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

P_VAL_THRESH = 1.0E-99
THRESH_MOD = 0.90

df_ivt = pd.read_csv(ivt_res_file, sep='\t')
df_ivt = df_ivt.rename(columns={'Unnamed: 0': 'site_index'})
df_wt = pd.read_csv(wt_res_file, sep='\t')
df_wt = df_wt.rename(columns={'Unnamed: 0': 'site_index'})

motifs = df_wt['ref_motif'].unique()

fig_hist = plt.figure(figsize=(16, 12))
axes_hist = fig_hist.subplots(2, 3).flatten()
fig_mod_ratio = plt.figure(figsize=(16, 12))
axes_mod_ratio = fig_mod_ratio.subplots(2, 3).flatten()

num_bins = 100
bin_edges = np.linspace(0, 1, num_bins + 1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

for subplot_ind, this_motif in enumerate(motifs):
    df_motif_ivt = df_ivt[
        (df_ivt['P_adjust'] <= P_VAL_THRESH)
        & (df_ivt['ref_motif']==this_motif)
        # & (df_ivt['pred_motif']==this_motif)
        ]
    df_motif_wt = df_wt[
        (df_wt['P_adjust'] <= P_VAL_THRESH)
        & (df_wt['ref_motif']==this_motif)
        # & (df_wt['pred_motif']==this_motif)
        ]

    ### histogram of mod prob per read ###
    bin_counts_ivt, _ = np.histogram(df_motif_ivt['mod_prob'], bins=bin_edges)
    norm_counts_ivt = bin_counts_ivt / np.sum(bin_counts_ivt)

    bin_counts_wt, _ = np.histogram(df_motif_wt['mod_prob'], bins=bin_edges)
    norm_counts_wt = bin_counts_wt / np.sum(bin_counts_wt)

    axes_hist[subplot_ind].plot(bin_centers, norm_counts_ivt, c='b', label='IVT, {} features'.format(len(df_motif_ivt)))
    axes_hist[subplot_ind].plot(bin_centers, norm_counts_wt, c='r', label='WT, {} features'.format(len(df_motif_wt)))

    axes_hist[subplot_ind].set_xlabel('Mod. Prob.', fontsize=15)
    axes_hist[subplot_ind].set_ylabel('Norm. frequency', fontsize=15)
    axes_hist[subplot_ind].set_title('{}'.format(this_motif), fontsize=20)
    axes_hist[subplot_ind].set_xlim([-0.05, 1.05])
    axes_hist[subplot_ind].set_ylim([0, 0.3])
    axes_hist[subplot_ind].legend(loc='upper left', fontsize=12)

    ### aggregate mod. ratio per site ###
    def calc_mod_ratio(in_df_motif, thresh_mod=0.5):
        df_motif_avg = pd.DataFrame()
        for site_index in in_df_motif['site_index'].unique():
            df_site = in_df_motif[in_df_motif['site_index']==site_index]
            num_features = len(df_site)
            mod_ratio = np.mean((df_site['mod_prob'].values)>=thresh_mod)
            new_row = df_site.iloc[0][:16]
            new_row['num_features'] = num_features
            new_row['mod_ratio'] = mod_ratio
            df_motif_avg = pd.concat([df_motif_avg, new_row.to_frame().T])
        df_motif_avg.reset_index(drop=True)
        return df_motif_avg

    df_motif_avg_ivt = calc_mod_ratio(df_motif_ivt, thresh_mod=THRESH_MOD)
    glori_ratio_ivt = np.float64(df_motif_avg_ivt['Ratio'])
    mod_ratio_ivt = np.float64(df_motif_avg_ivt['mod_ratio'])
    corr_ivt = np.corrcoef(glori_ratio_ivt, mod_ratio_ivt)[0, 1]

    df_motif_avg_wt = calc_mod_ratio(df_motif_wt, thresh_mod=THRESH_MOD)
    glori_ratio_wt = np.float64(df_motif_avg_wt['Ratio'])
    mod_ratio_wt = np.float64(df_motif_avg_wt['mod_ratio'])
    corr_wt = np.corrcoef(glori_ratio_wt, mod_ratio_wt)[0, 1]

    # axes_mod_ratio[subplot_ind].plot(glori_ratio_ivt, mod_ratio_ivt, 'bo', mfc='none', label='IVT, {} sites, corr. {:.2f}'.format(len(df_motif_avg_ivt), corr_ivt))
    # axes_mod_ratio[subplot_ind].plot(glori_ratio_wt, mod_ratio_wt, 'ro', mfc='none', label='WT, {} sites, corr. {:.2f}'.format(len(df_motif_avg_wt), corr_wt))
    axes_mod_ratio[subplot_ind].plot(glori_ratio_ivt, mod_ratio_ivt, 'b.', label='IVT, {} sites, corr. {:.2f}'.format(len(df_motif_avg_ivt), corr_ivt))
    axes_mod_ratio[subplot_ind].plot(glori_ratio_wt, mod_ratio_wt, 'r.', label='WT, {} sites, corr. {:.2f}'.format(len(df_motif_avg_wt), corr_wt))

    axes_mod_ratio[subplot_ind].plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'k--', alpha=0.5)
    axes_mod_ratio[subplot_ind].set_xlim([-0.05, 1.05])
    axes_mod_ratio[subplot_ind].set_ylim([-0.05, 1.05])
    axes_mod_ratio[subplot_ind].set_xlabel('GLORI mod. ratio', fontsize=15)
    axes_mod_ratio[subplot_ind].set_ylabel('ONT mod. ratio', fontsize=15)
    axes_mod_ratio[subplot_ind].set_title('{}'.format(this_motif), fontsize=15)
    axes_mod_ratio[subplot_ind].legend(loc='upper left', fontsize=12)

# plt.subplots_adjust(top=0.8)
fig_hist.tight_layout()
fig_hist.savefig(os.path.join(img_out, 'hist_modProbs_pValThresh{:.2E}_modThresh{:.2f}.png'.format(P_VAL_THRESH, THRESH_MOD)), bbox_inches='tight')
fig_mod_ratio.tight_layout()
fig_mod_ratio.savefig(os.path.join(img_out, 'corr_glori_modRatio_pValThresh{:.2E}_modThresh{:.2f}.png'.format(P_VAL_THRESH, THRESH_MOD)), bbox_inches='tight')
# plt.close('all')