import os
HOME = os.path.expanduser('~')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_multiple_motifs_minCov10_plusMinusReads.tsv'
img_out = os.path.join(HOME, 'img_out', 'MAFIA_MaxAbs_multiple')
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

P_VAL_THRESH = 1.0E-10
NUM_READS_THRESH = 20
df_all = pd.read_csv(df_file, sep='\t')
df = df_all[(df_all['P_adjust']<=P_VAL_THRESH) * (df_all['num_test_features']>=NUM_READS_THRESH)]
# df = df_all

motifs = df['motif'].unique()

for this_motif in motifs:
    df_motif = df[df['motif']==this_motif]
    corr = np.corrcoef(df_motif['Ratio'], df_motif['mod_ratio'])[0, 1]

    plt.figure(figsize=(6, 6))
    plt.plot(df_motif['Ratio'], df_motif['mod_ratio'], 'o', mfc='none')
    # plt.plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'r-', alpha=0.5)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('GLORI (p_val_thresh {:.2E})'.format(P_VAL_THRESH), fontsize=15)
    plt.ylabel('Predicted mod. ratio', fontsize=15)
    plt.title('{}, {} sites >= {} reads\nCorrelation {:.2f}'.format(this_motif, df_motif.shape[0], NUM_READS_THRESH, corr), fontsize=20)
    plt.show()

    plt.savefig(os.path.join(img_out, 'corr_glori_modRatio_{}_minNumReads{}_pValThresh{:.2E}.png'.format(this_motif, NUM_READS_THRESH, P_VAL_THRESH)), bbox_inches='tight')
    plt.close()