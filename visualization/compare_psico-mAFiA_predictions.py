import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import os

THRESH_CONF = 80
THRESH_COV = 50
pred_ds = 'Mettl3-KO_rep1'
bid_seq_calibrated = False
restrict_motifs = True
# comp_ds = 'BID-Seq'
# mod_type = 'psi'
comp_ds = 'WT_P2'
mod_type = 'm6A'

dict_ds = {
    'RNA004_HEK293_WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/RNA004/dorado/RNA004_HEK293_WT_RTA.mAFiA.bed',

    'WT_P2': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/chrALL.mAFiA.sites.bed',

    '100_WT_0_IVT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HEK293/100_WT_0_IVT/chrALL.mAFiA.sites.bed',
    '75_WT_25_IVT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HEK293/75_WT_25_IVT/chrALL.mAFiA.sites.bed',
    '50_WT_50_IVT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HEK293/50_WT_50_IVT/chrALL.mAFiA.sites.bed',
    '25_WT_75_IVT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HEK293/25_WT_75_IVT/chrALL.mAFiA.sites.bed',
    '0_WT_100_IVT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HEK293/0_WT_100_IVT/chrALL.mAFiA.sites.bed',
    'Mettl3-KO_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HEK293T-Mettl3-KO/rep1//chrALL.mAFiA.sites.bed',

    'HeLa_WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/HeLa/chrALL.mAFiA.sites.bed',
    # 'HeLa_SRR28796313': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HeLa_SRR28796313/chrALL.mAFiA.sites.bed',

    'siCtrl_input_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siCtrl_input_rep1/chrALL.mAFiA.sites.bed',
    'siMETTL3_input_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siMETTL3_input_rep1/chrALL.mAFiA.sites.bed',
    'siTRUB1_input_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siTRUB1_input_rep1/chrALL.mAFiA.sites.bed',
    'siCtrl_input_merged': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siCtrl_input_rep1/chrALL.mAFiA.sites.bed',
    'siMETTL3_input_merged': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siMETTL3_input_rep1/chrALL.mAFiA.sites.bed',
    'siTRUB1_input_merged': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siTRUB1_input_rep1/chrALL.mAFiA.sites.bed',

    'TAC_ctrl': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/TAC/40-34/chrALL.mAFiA.sites.bed',
    'SHAM_day1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/TAC/40-31/chrALL.mAFiA.sites.bed',
    'A1_WT_CD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/A1_WT_CD/chrALL.mAFiA.sites.bed',
    'B1_M3KO_CD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/B1_M3KO_CD/chrALL.mAFiA.sites.bed',
    'C1_WT_WD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/C1_WT_WD/chrALL.mAFiA.sites.bed',
    'D1_M3KO_WD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/D1_M3KO_WD/chrALL.mAFiA.sites.bed',

    'GLORI': '/home/adrian/Data/GLORI/bed_files/GLORI.chrALL.tsv',
    'BID-Seq': '/home/adrian/Data/BID_seq/BID_seq_HEK293T.bed',
    'BID-Seq_mouse_heart': '/home/adrian/Data/BID_seq/BID_seq_mouse_heart.bed',
    'PRAISE': '/home/adrian/Data/PRAISE/PRAISE_HEK293T_span1-3.bed',
    'BACS': '/home/adrian/Data/BACS/BACS_HeLa_WT.bed'
}

df_pred = pd.read_csv(
    dict_ds[pred_ds],
    sep='\t',
    dtype={'chrom': str}
)
if 'confidence' in df_pred.keys():
    df_pred = df_pred[
        (df_pred['confidence']>=THRESH_CONF)
        * (df_pred['coverage']>=THRESH_COV)
    ]
else:
    df_pred = df_pred[df_pred['coverage']>=THRESH_COV]

df_comp = pd.read_csv(
    dict_ds[comp_ds],
    sep='\t',
    dtype={'chrom': str}
)
if 'confidence' in df_comp.keys():
    df_comp = df_comp[df_comp['confidence']>=THRESH_CONF]

if bid_seq_calibrated:
    df_comp.rename(columns={'BID-Seq': 'modRatio'}, inplace=True)
elif comp_ds in ['GLORI', 'BID-Seq', 'BID-Seq_mouse_heart', 'PRAISE', 'BACS']:
    df_comp.rename(columns={'score': 'modRatio'}, inplace=True)

img_out = f'/home/adrian/img_out/psi-co-mAFiA/{pred_ds}'
os.makedirs(img_out, exist_ok=True)

def scatter_plot(df_in, key_x, key_y, mod_type, fig_name, corr=None):
    plt.figure(figsize=(5, 5))
    plt.plot(df_in[key_x], df_in[key_y], '.', alpha=0.5)
    plt.plot([0, 100], [0, 100], linestyle='--', c='r', alpha=0.5)
    plt.xlim([-1, 101])
    plt.ylim([-1, 101])
    plt.xticks(np.linspace(0, 100, 5))
    plt.yticks(np.linspace(0, 100, 5))
    plt.xlabel("$S_{{{}}}$".format('-'.join(key_x.split('_')[1:])), fontsize=10)
    plt.ylabel("$S_{{{}}}$".format('-'.join(key_y.split('_')[1:])), fontsize=10)
    if corr is not None:
        title = f'{len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%, corr. {corr:.2f}'
    else:
        title = f'{len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%'
    plt.suptitle(title, fontsize=12)
    plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')

def scatter_plot_by_motif(df_in, key_x, key_y, mod_type, ordered_motifs, num_row, num_col, fig_name, thresh_err=25, calc_error=False):
    motif_counts = Counter(df_in['ref5mer'])
    if calc_error:
        motif_err_rate = {
            motif: (df_in[df_in['ref5mer'] == motif][key_y] >= thresh_err).sum() / motif_counts[motif] for
            motif in motif_counts.keys()}

    plt.figure(figsize=(num_col*2, num_row*2))
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
            plt.xlabel("$S_{{{}}}$".format('-'.join(key_x.split('_')[1:])), fontsize=10)
        if ind % num_col == 0:
            plt.ylabel("$S_{{{}}}$".format('-'.join(key_y.split('_')[1:])), fontsize=10)
    plt.suptitle(f'{len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%', fontsize=12)
    plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')

if mod_type=='m6A':
    motifs = [
        'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
        'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
         'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
    ]
    num_rows = 3
    num_cols = 6
elif mod_type=='psi':
    motifs = [
        'GTTCA', 'GTTCG', 'GTTCC', 'GTTCT', 'TGTAG', 'TGTGG',
        'GGTGG', 'ATTTG', 'GATGC', 'AGTGG', 'CATAA', 'CATCC',
        'CCTCC', 'CTTTA', 'GGTCC', 'TATAA'
    ]
    num_rows = 3
    num_cols = 6

def calc_correlation(in_df):
    in_array = in_df[[f'modRatio_{pred_ds}', f'modRatio_{comp_ds}']].values
    out_num_sites = len(in_array)
    out_corr = np.corrcoef(in_array.T)[0, 1]
    return round(out_corr, 3), out_num_sites

### compare to Bid-Seq ###
df_comp_pred = pd.merge(df_comp, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=[f'_{comp_ds}', f'_{pred_ds}'])
df_comp_pred = df_comp_pred[df_comp_pred['name']==mod_type]
if restrict_motifs:
    df_comp_pred_sel = df_comp_pred[df_comp_pred['ref5mer'].isin(motifs)]
else:
    df_comp_pred_sel = df_comp_pred
corr, num_sites = calc_correlation(df_comp_pred_sel)
with open(os.path.join(img_out, f'corr_{mod_type}_pred_vs_{comp_ds}_conf{THRESH_CONF}_cov{THRESH_COV}.txt'), 'w') as f_out:
    f_out.write('num_sites' + '\t' + 'correlation' + '\n')
    f_out.write(str(num_sites) + '\t' + str(corr) + '\t' + '\n')
scatter_plot_by_motif(df_comp_pred_sel, f'modRatio_{comp_ds}', f'modRatio_{pred_ds}', mod_type, motifs, num_rows, num_cols, f'{mod_type}_pred_vs_{comp_ds}_conf{THRESH_CONF}_cov{THRESH_COV}.png')
if restrict_motifs:
    out_filename = f'{mod_type}_pred_vs_{comp_ds}_combined_conf{THRESH_CONF}_cov{THRESH_COV}_restrict_motifs.png'
else:
    out_filename = f'{mod_type}_pred_vs_{comp_ds}_combined_conf{THRESH_CONF}_cov{THRESH_COV}.png'
scatter_plot(df_comp_pred_sel, f'modRatio_{comp_ds}', f'modRatio_{pred_ds}', mod_type, out_filename, corr=corr)
