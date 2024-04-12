import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import os

chrom = 'ALL'
THRESH_CONF = 80

# ds = '100_WT_0_IVT'
ds = 'HeLa'
# ds = '40-34'
# ds = 'A2_WT_CD'

# comp = 'BID-Seq_observed'
comp = 'BID-Seq_calibrated'
# comp = 'PRAISE'
# comp = 'Psi-Seq'

########################################################################################################################
if ds == '100_WT_0_IVT':
    pred_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/{ds}/chrALL.mAFiA.sites.bed'
    if comp in ['BID-Seq_observed', 'BID-Seq_calibrated']:
        comp_file = '/home/adrian/Data/BID_seq/BID_seq_HEK293T.bed'
    elif comp == 'PRAISE':
        comp_file = '/home/adrian/Data/PRAISE/PRAISE_HEK293T.bed'
    elif comp == 'Psi-Seq':
        comp_file = '/home/adrian/Data/Psi-Seq/Psi-Seq_HEK293.bed'
elif ds == 'HeLa':
    pred_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/{ds}/chrALL.mAFiA.sites.bed'
    comp_file = '/home/adrian/Data/BID_seq/BID_seq_HeLa.bed'
elif ds == '40-34':
    pred_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/mouse_heart/Dewenter_TAC_Backs_lab/{ds}/chrALL.mAFiA.sites.bed'
    comp_file = '/home/adrian/Data/BID_seq/BID_seq_mouse_heart.bed'
elif ds in ['A1_WT_CD', 'A2_WT_CD']:
    pred_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/mouse_heart/Federica_Accornero/{ds}/chrALL.mAFiA.sites.bed'
    comp_file = '/home/adrian/Data/BID_seq/BID_seq_mouse_heart.bed'

########################################################################################################################

df_pred = pd.read_csv(
    pred_file,
    sep='\t',
    dtype={'chrom': str}
)

df_comp = pd.read_csv(
    comp_file,
    sep='\t',
    dtype={'chrom': str}
)
if comp in ['PRAISE', 'Psi-Seq']:
    df_comp.rename(columns={'score': 'modRatio'}, inplace=True)
elif comp=='BID-Seq_observed':
    df_comp.rename(columns={'score': 'modRatio'}, inplace=True)
elif comp=='BID-Seq_calibrated':
    df_comp.rename(columns={'BID-Seq': 'modRatio'}, inplace=True)

img_out = f'/home/adrian/img_out/psi-co-mAFiA/{ds}'
os.makedirs(img_out, exist_ok=True)

def scatter_plot(df_in, key_x, key_y, mod_type, fig_name, corr=None):
    plt.figure(figsize=(5, 5))
    plt.plot(df_in[key_x], df_in[key_y], 'o', alpha=0.5)
    plt.plot([0, 100], [0, 100], linestyle='--', c='r', alpha=0.5)
    plt.xlim([-1, 101])
    plt.ylim([-1, 101])
    plt.xticks(np.linspace(0, 100, 5))
    plt.yticks(np.linspace(0, 100, 5))
    plt.xlabel("$S_{{{}}}$".format(key_x.split('_')[1]), fontsize=10)
    plt.ylabel("$S_{{{}}}$".format((key_y.split('_')[1])), fontsize=10)
    if corr is not None:
        title = f'{ds} chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%, {comp}\ncorr. {corr:.2f}'
    else:
        title = f'{ds} chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%, {comp}'
    plt.suptitle(title, fontsize=12)
    plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')

def scatter_plot_by_motif(df_in, key_x, key_y, mod_type, ordered_motifs, num_row, num_col, fig_name, thresh_err=25, calc_error=False):
    motif_counts = Counter(df_in['ref5mer'])
    if calc_error:
        motif_err_rate = {
            motif: (df_in[df_in['ref5mer'] == motif][key_y] >= thresh_err).sum() / motif_counts[motif] for
            motif in motif_counts.keys()}

    plt.figure(figsize=(12, 6))
    for ind, motif in enumerate(ordered_motifs):
        # if motif not in df_in['ref5mer'].values:
        #     continue
        sub_df = df_in[df_in['ref5mer'] == motif]
        plt.subplot(num_row, num_col, ind + 1)
        # plt.scatter(sub_df[key_x], sub_df[key_y], s=2, alpha=0.5)
        plt.plot(sub_df[key_x], sub_df[key_y], '.', alpha=0.5)
        if calc_error:
            plt.text(1, 90, f"{motif.replace('T', 'U')} ({motif_err_rate[motif]:.2f})", fontsize=10)
            plt.axhline(y=thresh_err, linestyle='--', c='r')
        else:
            plt.text(1, 90, f"{motif.replace('T', 'U')}", fontsize=10)
        plt.plot([0, 100], [0, 100], linestyle='--', c='r', alpha=0.5)
        plt.xlim([-1, 101])
        plt.ylim([-1, 101])
        if ind >= (num_row-1)*num_col:
            plt.xticks(np.linspace(0, 100, 5))
        else:
            plt.xticks(np.linspace(0, 100, 5), [])
        if ind % num_col == 0:
            plt.yticks(np.linspace(0, 100, 5))
        else:
            plt.yticks(np.linspace(0, 100, 5), [])
        if ind >= ((num_row - 1) * num_col):
            plt.xlabel("$S_{{{}}}$".format(key_x.split('_')[1]), fontsize=10)
        if ind % num_col == 0:
            plt.ylabel("$S_{{{}}}$".format((key_y.split('_')[1])), fontsize=10)
    plt.suptitle(f'{ds} chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%, {comp}', fontsize=12)
    plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')

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

# m6A_motifs = [
#     'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
#     'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
#      'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
# ]

def calc_correlation(in_df):
    in_array = in_df[['modRatio_pred', f'modRatio_{comp}']].values
    out_num_sites = len(in_array)
    out_corr = np.corrcoef(in_array.T)[0, 1]
    return round(out_corr, 3), out_num_sites

### compare to Bid-Seq ###
df_comp_pred = pd.merge(df_comp, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=[f'_{comp}', '_pred'])
df_comp_pred = df_comp_pred[df_comp_pred['name']=='psi']
df_comp_pred_sel = df_comp_pred[
    (df_comp_pred['confidence']>=THRESH_CONF)
]
corr, num_sites = calc_correlation(df_comp_pred_sel)
with open(os.path.join(img_out, f'corr_psi_pred_vs_{comp}_conf{THRESH_CONF}.txt'), 'w') as f_out:
    f_out.write('num_sites' + '\t' + 'correlation' + '\n')
    f_out.write(str(num_sites) + '\t' + str(corr) + '\t' + '\n')
scatter_plot_by_motif(df_comp_pred_sel, f'modRatio_{comp}', 'modRatio_pred', 'psi', psi_motifs, 3, 5, f'psi_pred_vs_{comp}_conf{THRESH_CONF}.png')
scatter_plot(df_comp_pred_sel, f'modRatio_{comp}', 'modRatio_pred', 'psi', f'psi_pred_vs_{comp}_combined_conf{THRESH_CONF}.png', corr=corr)
