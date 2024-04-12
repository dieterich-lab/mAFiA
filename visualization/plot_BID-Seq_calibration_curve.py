import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import numpy as np

bid_file = '/home/adrian/Data/BID_seq/41587_2022_1505_MOESM3_ESM.xlsx'
df_in = pd.read_excel(bid_file, sheet_name='Supplementary Table 1', skiprows=range(3))

img_out = '/home/adrian/img_out/BID-Seq'
os.makedirs(img_out, exist_ok=True)

psi_motifs = [
    'AGTGG',
    'ATTTG',
    'CATAA',
    'CATCC',
    'CCTCC',
    'GATGC',
    'GGTCC',
    'GGTGG',
    'GTTCA',
    'GTTCC',
    'GTTCG',
    'GTTCT',
    'TATAA',
    'TGTAG',
    'TGTGG'
]

def calibrate(in_y, in_pars):
    return (in_y - in_pars[1]) / (in_pars[2] - in_pars[0]*in_pars[2] - pars[1] + pars[0]*in_y) * 100.0


num_row = 3
num_col = 5
vec_y = np.linspace(0, 1, 101)
plt.figure(figsize=(12, 6))
for ind, this_motif in enumerate(psi_motifs):
    pars = df_in[df_in['Motif']==this_motif][['fit_A', 'fit_B', 'fit_R']].values[0]
    vec_x = calibrate(vec_y, pars)

    plt.subplot(3, 5, ind+1)
    plt.plot(vec_y, vec_x, label=this_motif)
    plt.legend(loc='lower right', fontsize=10)

    if ind >= (num_row - 1) * num_col:
        plt.xticks(np.linspace(0, 1, 5))
    else:
        plt.xticks(np.linspace(0, 1, 5), [])
    # if ind % num_col == 0:
    #     plt.yticks(np.linspace(0, 100, 5))
    # else:
    #     plt.yticks(np.linspace(0, 100, 5), [])

    if ind >= ((num_row - 1) * num_col):
        plt.xlabel('Observed deletion frac.', fontsize=10)
    if ind % num_col == 0:
        plt.ylabel('Calibrated %', fontsize=10)
plt.savefig(os.path.join(img_out, f'BID-Seq_calibration_curve_{len(psi_motifs)}.png'), bbox_inches='tight')