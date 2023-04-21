import os
HOME = os.path.expanduser('~')
import argparse
import pandas as pd
import numpy as np
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# parser = argparse.ArgumentParser()
# parser.add_argument('--df_file')
# args = parser.parse_args()
# df_file = args.df_file
df_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_random_ligation_A_m6A_modProbThresh0.5.tsv'
img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename(df_file).rstrip('.tsv'))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

P_VAL_THRESH = 1.0E-99
# NUM_READS_THRESH = 20
df_all = pd.read_csv(df_file, sep='\t')
df_all = df_all.rename(columns={'Unnamed: 0': 'site_index'})
# df = df_all

motifs = df_all['ref_motif'].unique()

fig_hist = plt.figure(figsize=(16, 12))
axes_hist = fig_hist.subplots(2, 3).flatten()
fig_mod_ratio = plt.figure(figsize=(16, 12))
axes_mod_ratio = fig_mod_ratio.subplots(2, 3).flatten()

for subplot_ind, this_motif in enumerate(motifs):
    df = df_all[df_all['P_adjust'] <= P_VAL_THRESH]
    df_motif = df[df['ref_motif']==this_motif]

    ### histogram of mod prob ###
    axes_hist[subplot_ind].hist(df_motif['mod_prob'], bins=50, range=[0, 1])
    axes_hist[subplot_ind].set_xlabel('Mod. Prob.', fontsize=15)
    axes_hist[subplot_ind].set_ylabel('Counts', fontsize=15)
    axes_hist[subplot_ind].set_title('{}, {} features'.format(this_motif, len(df_motif)), fontsize=20)

    ### average per site ###
    df_motif_avg = pd.DataFrame()
    for site_index in df_motif['site_index'].unique():
        df_site = df_motif[df_motif['site_index']==site_index]
        num_features = len(df_site)
        mod_ratio = df_site['mod_prob'].mean()
        new_row = df_site.iloc[0][:16]
        new_row['num_features'] = num_features
        new_row['mod_ratio'] = mod_ratio
        df_motif_avg = pd.concat([df_motif_avg, new_row.to_frame().T])
    df_motif_avg.reset_index(drop=True)

    glori_ratio = np.float64(df_motif_avg['Ratio'])
    mod_ratio = np.float64(df_motif_avg['mod_ratio'])
    corr = np.corrcoef(glori_ratio, mod_ratio)[0, 1]

    axes_mod_ratio[subplot_ind].plot(glori_ratio, mod_ratio, 'o', mfc='none')
    axes_mod_ratio[subplot_ind].plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'k--', alpha=0.5)
    axes_mod_ratio[subplot_ind].set_xlim([-0.05, 1.05])
    axes_mod_ratio[subplot_ind].set_ylim([-0.05, 1.05])
    axes_mod_ratio[subplot_ind].set_xlabel('GLORI (p-val$\leq${:.2E})'.format(P_VAL_THRESH), fontsize=15)
    axes_mod_ratio[subplot_ind].set_ylabel('ONT mod. ratio', fontsize=15)
    axes_mod_ratio[subplot_ind].set_title('{}, {} sites'.format(this_motif, df_motif_avg.shape[0]) + '\nCorrelation {:.2f}'.format(corr), fontsize=15)
# plt.subplots_adjust(top=0.8)
fig_hist.tight_layout()
fig_hist.savefig(os.path.join(img_out, 'hist_modProbs_pValThresh{:.2E}.png'.format(P_VAL_THRESH)), bbox_inches='tight')
fig_mod_ratio.tight_layout()
fig_mod_ratio.savefig(os.path.join(img_out, 'corr_glori_modRatio_pValThresh{:.2E}.png'.format(P_VAL_THRESH)), bbox_inches='tight')
plt.close('all')