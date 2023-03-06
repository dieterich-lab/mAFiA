import os
HOME = os.path.expanduser('~')
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--df_file')
args = parser.parse_args()
df_file = args.df_file
# df_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_multipleMotifs_noNorm_g3.tsv'
img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename(df_file).rstrip('.tsv'))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

P_VAL_THRESH = 1.0E-99
# NUM_READS_THRESH = 20
df_all = pd.read_csv(df_file, sep='\t')
# df = df_all

motifs = df_all['motif'].unique()

for this_motif in motifs:
    plt.figure(figsize=(18, 6))
    for subplot_ind, coverage_thresh in enumerate([10, 20, 50]):
        df = df_all[(df_all['P_adjust'] <= P_VAL_THRESH) * (df_all['num_test_features'] >= coverage_thresh)]
        df_motif = df[df['motif']==this_motif]
        corr = np.corrcoef(df_motif['Ratio'], df_motif['mod_ratio'])[0, 1]

        plt.subplot(1, 3, subplot_ind+1)
        plt.plot(df_motif['Ratio'], df_motif['mod_ratio'], 'o', mfc='none')
        # plt.plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'r-', alpha=0.5)
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('GLORI (p-val$\leq${:.2E})'.format(P_VAL_THRESH), fontsize=15)
        if subplot_ind==0:
            plt.ylabel('ONT mod. ratio', fontsize=15)
        plt.title(r'{} sites with ONT coverage$\geq${}'.format(df_motif.shape[0], coverage_thresh) + '\nCorrelation {:.2f}'.format(corr), fontsize=15)

    plt.subplots_adjust(top=0.8)
    plt.suptitle(this_motif, fontsize=20)

    plt.savefig(os.path.join(img_out, 'corr_glori_modRatio_{}_pValThresh{:.2E}.png'.format(this_motif, P_VAL_THRESH)), bbox_inches='tight')
    plt.close()