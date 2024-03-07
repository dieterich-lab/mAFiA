import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import Counter
import pysam
import numpy as np
from tqdm import tqdm
import os
from Bio.Seq import Seq

chrom = 'ALL'
ds = '40-34'

THRESH_CONF = 80

pred_file = f'/home/adrian/Data/Dewenter_TAC_Backs_lab/achan/psico-mAFiA_results/{ds}/chrALL.mAFiA.sites.bed'
bid_file = '/home/adrian/Data/BID_seq/BID_seq_mouse_heart.bed'
# bid_file = '/home/adrian/Data/BID_seq/BID_seq_HEK293T.bed'

df_pred = pd.read_csv(
    pred_file,
    sep='\t',
    dtype={'chrom': str}
)

df_bid = pd.read_csv(
    bid_file,
    sep='\t',
    dtype={'chrom': str}
)
df_bid.rename(columns={'score': 'modRatio'}, inplace=True)

img_out = '/home/adrian/img_out/mouse_heart'
os.makedirs(img_out, exist_ok=True)

def scatter_plot(df_in, key_x, key_y, mod_type, fig_name):
    plt.figure(figsize=(5, 5))
    plt.plot(df_in[key_x], df_in[key_y], 'o', alpha=0.5)
    plt.xlim([-1, 101])
    plt.ylim([-1, 101])
    plt.xticks(np.linspace(0, 100, 5))
    plt.yticks(np.linspace(0, 100, 5))
    plt.xlabel("$S_{{{}}}$".format((key_x.split('_')[1]).upper()), fontsize=10)
    plt.ylabel("$S_{{{}}}$".format((key_y.split('_')[1]).upper()), fontsize=10)
    plt.suptitle(f'{ds} chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%', fontsize=12)
    plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')

def scatter_plot_by_motif(df_in, key_x, key_y, mod_type, ordered_motifs, num_row, num_col, fig_name, thresh_err=25, calc_error=False):
    motif_counts = Counter(df_in['ref5mer'])
    if calc_error:
        motif_err_rate = {
            motif: (df_in[df_in['ref5mer'] == motif][key_y] >= thresh_err).sum() / motif_counts[motif] for
            motif in motif_counts.keys()}

    plt.figure(figsize=(12, 6))
    for ind, motif in enumerate(ordered_motifs):
        if motif not in df_in['ref5mer'].values:
            continue
        sub_df = df_in[df_in['ref5mer'] == motif]
        plt.subplot(num_row, num_col, ind + 1)
        # plt.scatter(sub_df[key_x], sub_df[key_y], s=2, alpha=0.5)
        plt.plot(sub_df[key_x], sub_df[key_y], '.', alpha=0.5)
        if calc_error:
            plt.text(1, 90, f"{motif.replace('T', 'U')} ({motif_err_rate[motif]:.2f})", fontsize=10)
            plt.axhline(y=thresh_err, linestyle='--', c='r')
        else:
            plt.text(1, 90, f"{motif.replace('T', 'U')}", fontsize=10)
        plt.xlim([-1, 101])
        plt.ylim([-1, 101])
        plt.xticks(np.linspace(0, 100, 5))
        plt.yticks(np.linspace(0, 100, 5))
        if ind >= ((num_row - 1) * num_col):
            plt.xlabel("$S_{{{}}}$".format((key_x.split('_')[1]).upper()), fontsize=10)
        if ind % num_col == 0:
            plt.ylabel("$S_{{{}}}$".format((key_y.split('_')[1]).upper()), fontsize=10)
    plt.suptitle(f'{ds} chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%', fontsize=12)
    plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')

psi_motifs = [
    'AGTGG',
    'GGTCC',
    'GGTGG',
    'GTTCA',
    'GTTCC',
    'GTTCG',
    'TGTAG',
    'TGTGG'
]

m6A_motifs = [
    'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
    'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
     'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
]

### compare to Bid-Seq ###
df_bid_pred = pd.merge(df_bid, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=['_bid', '_pred'])
df_bid_pred_sel = df_bid_pred[
    (df_bid_pred['confidence']>=THRESH_CONF)
]
scatter_plot_by_motif(df_bid_pred_sel, 'modRatio_bid', 'modRatio_pred', 'psi', psi_motifs, 2, 4, f'{ds}_psi_pred_vs_BID_conf{THRESH_CONF}.png')
scatter_plot(df_bid_pred_sel, 'modRatio_bid', 'modRatio_pred', 'psi', f'{ds}_psi_pred_vs_BID_combined_conf{THRESH_CONF}.png')
