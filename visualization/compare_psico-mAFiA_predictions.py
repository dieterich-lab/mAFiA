import pandas as pd
from collections import Counter
import numpy as np
import os
from tqdm import tqdm
import matplotlib as mpl
#######################################################################
# cm = 1/2.54  # centimeters in inches
# gr = 1.618
# # mpl.rcParams['figure.dpi'] = 600
# # mpl.rcParams['savefig.dpi'] = 600
# mpl.rcParams['font.size'] = 5
# mpl.rcParams['legend.fontsize'] = 5
# mpl.rcParams['xtick.labelsize'] = 5
# mpl.rcParams['ytick.labelsize'] = 5
# mpl.rcParams['xtick.major.size'] = 1.5
# mpl.rcParams['ytick.major.size'] = 1.5
# mpl.rcParams['lines.linewidth'] = 0.5
# mpl.rcParams['font.family'] = 'Arial'
# FMT = 'pdf'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################
import matplotlib.pyplot as plt

THRESH_CONF = 80
THRESH_COV = 20
pred_ds = 'CM_M3KO_rep2'
bid_seq_calibrated = False
restrict_motifs = None
comp_ds = 'CM_WT_rep2'
mod_type = 'm6A'

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

dict_ds = {
    'RNA004_HEK293_WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/RNA004/dorado/RNA004_HEK293_WT_RTA.mAFiA.bed',

    'HEK293_WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/chrALL.mAFiA.sites.bed',

    '100WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/100WT/chrALL.mAFiA.sites.bed',
    '75WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/75WT/chrALL.mAFiA.sites.bed',
    '50WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/50WT/chrALL.mAFiA.sites.bed',
    '25WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/25WT/chrALL.mAFiA.sites.bed',
    '0WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/0WT/chrALL.mAFiA.sites.bed',

    'METTL3_KO_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293T_Mettl3_KO/rep1/chrALL.mAFiA.sites.bed',
    'METTL3_KO_rep2': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293T_Mettl3_KO/rep2/chrALL.mAFiA.sites.bed',
    'METTL3_KO_rep3': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293T_Mettl3_KO/rep3/chrALL.mAFiA.sites.bed',
    'METTL3_KO': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293T_Mettl3_KO/merged/chrALL.mAFiA.sites.bed',

    'TRUB1_OE_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293_TRUB1_OE/rep1/chrALL.mAFiA.sites.bed',
    'TRUB1_OE_rep2': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293_TRUB1_OE/rep2/chrALL.mAFiA.sites.bed',
    'TRUB1_OE': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293_TRUB1_OE/merged/chrALL.mAFiA.sites.bed',

    'HeLa_NanoSPA_WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/NanoSPA/HeLa_WT/chrALL.mAFiA.sites.bed',
    'HeLa_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HeLa/rep1/chrALL.mAFiA.sites.bed',
    'HeLa_rep2': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HeLa/rep2/chrALL.mAFiA.sites.bed',
    'HeLa_rep3': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HeLa/rep3/chrALL.mAFiA.sites.bed',
    # 'HeLa_SRR28796313': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HeLa_SRR28796313/chrALL.mAFiA.sites.bed',

    'siCtrl_input_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siCtrl_input_rep1/chrALL.mAFiA.sites.bed',
    'siMETTL3_input_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siMETTL3_input_rep1/chrALL.mAFiA.sites.bed',
    'siTRUB1_input_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siTRUB1_input_rep1/chrALL.mAFiA.sites.bed',
    'siCtrl_input_merged': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/NanoSPA/HEK_siCtrl_input_rep1/chrALL.mAFiA.sites.bed',
    'METTL3_KD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/NanoSPA/HEK_siMETTL3_input_merged//chrALL.mAFiA.sites.bed',
    'TRUB1_KD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/NanoSPA/HEK_siTRUB1_input_merged/chrALL.mAFiA.sites.bed',

    'TAC_ctrl': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/TAC/40-34/chrALL.mAFiA.sites.bed',
    'SHAM_day1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/TAC/40-31/chrALL.mAFiA.sites.bed',
    'A1_WT_CD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/A1_WT_CD/chrALL.mAFiA.sites.bed',
    'B1_M3KO_CD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/B1_M3KO_CD/chrALL.mAFiA.sites.bed',
    'C1_WT_WD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/C1_WT_WD/chrALL.mAFiA.sites.bed',
    'D1_M3KO_WD': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/Diet/D1_M3KO_WD/chrALL.mAFiA.sites.bed',

    'CM_M3KO_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/CM/M3KO_rep1/chrALL.mAFiA.sites.bed',
    'CM_WT_rep1': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/CM/WT_rep1/chrALL.mAFiA.sites.bed',
    'CM_M3KO_rep2': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/CM/M3KO_rep2/chrALL.mAFiA.sites.bed',
    'CM_WT_rep2': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/CM/WT_rep2/chrALL.mAFiA.sites.bed',

    'CM_M3KO': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/CM/M3KO_merged/chrALL.mAFiA.sites.bed',
    'CM_WT': '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/CM/WT_merged/chrALL.mAFiA.sites.bed',

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
        if restrict_motifs:
            title = f'{restrict_motifs}\n{len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%, corr. {corr:.2f}'
        else:
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
    if restrict_motifs=='GUUCN':
        motifs = [
            'GTTCA', 'GTTCG', 'GTTCC', 'GTTCT'
        ]
    elif restrict_motifs=='Others':
        motifs = [
        'TGTAG', 'TGTGG', 'GGTGG', 'ATTTG',
        'GATGC', 'AGTGG', 'CATAA', 'CATCC',
        'CCTCC', 'CTTTA', 'GGTCC', 'TATAA'
    ]
    else:
        motifs = [
            'GTTCA', 'GTTCG', 'GTTCC', 'GTTCT',
            'TGTAG', 'TGTGG', 'GGTGG', 'ATTTG',
            'GATGC', 'AGTGG', 'CATAA', 'CATCC',
            'CCTCC', 'CTTTA', 'GGTCC', 'TATAA'
        ]
    num_rows = 4
    num_cols = 4

def calc_correlation(in_df):
    in_array = in_df[[f'modRatio_{pred_ds}', f'modRatio_{comp_ds}']].values
    out_num_sites = len(in_array)
    out_corr = np.corrcoef(in_array.T)[0, 1]
    return round(out_corr, 3), out_num_sites

### compare to Bid-Seq ###
df_comp_pred = pd.merge(df_comp, df_pred, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'], suffixes=[f'_{comp_ds}', f'_{pred_ds}'])
df_comp_pred_sel = df_comp_pred[df_comp_pred['name']==mod_type]
if restrict_motifs:
    df_comp_pred_sel = df_comp_pred_sel[df_comp_pred_sel['ref5mer'].isin(motifs)]
else:
    df_comp_pred_sel = df_comp_pred_sel
corr, num_sites = calc_correlation(df_comp_pred_sel)
with open(os.path.join(img_out, f'corr_{mod_type}_pred_vs_{comp_ds}_conf{THRESH_CONF}_cov{THRESH_COV}.txt'), 'w') as f_out:
    f_out.write('num_sites' + '\t' + 'correlation' + '\n')
    f_out.write(str(num_sites) + '\t' + str(corr) + '\t' + '\n')
scatter_plot_by_motif(df_comp_pred_sel, f'modRatio_{comp_ds}', f'modRatio_{pred_ds}', mod_type, motifs, num_rows, num_cols, f'{mod_type}_pred_vs_{comp_ds}_conf{THRESH_CONF}_cov{THRESH_COV}.png')
if restrict_motifs:
    out_filename = f'{mod_type}_pred_vs_{comp_ds}_combined_conf{THRESH_CONF}_cov{THRESH_COV}_restrict_motifs_{restrict_motifs}.png'
else:
    out_filename = f'{mod_type}_pred_vs_{comp_ds}_combined_conf{THRESH_CONF}_cov{THRESH_COV}.png'
scatter_plot(df_comp_pred_sel, f'modRatio_{comp_ds}', f'modRatio_{pred_ds}', mod_type, out_filename, corr=corr)

### histogram of deltaS ###
bin_max = 100
plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    sub_df = df_comp_pred[df_comp_pred['name']==this_mod]
    delta = sub_df[f'modRatio_{pred_ds}'] - sub_df[f'modRatio_{comp_ds}']
    plt.hist(delta[delta>=0], bins=bin_max, range=[0, bin_max], histtype='step', facecolor='b', label=f'$\Delta S \geq 0$')
    plt.hist(-delta[delta<0], bins=bin_max, range=[0, bin_max], histtype='step', facecolor='r', label=f'$\Delta S < 0$')
    plt.legend(loc='upper right')
    # plt.axvline(x=0, c='r')
    plt.xlabel(f'$Abs(\Delta S_{{{dict_mod_display[this_mod]}}})$', fontsize=12)
    plt.ylabel('Site Counts', fontsize=12)
    plt.yscale('log')
plt.suptitle(f'{pred_ds} cf. {comp_ds}', fontsize=15)
plt.savefig(os.path.join(img_out, f'deltaS_{pred_ds}_cf_{comp_ds}.png'), bbox_inches='tight')

### correlation vs distance ###
df_gene_bed = pd.read_csv('/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/gene.GRCh38.102.bed', sep='\t')

df_comp_pred['delta'] = df_comp_pred[f'modRatio_{pred_ds}'] - df_comp_pred[f'modRatio_{comp_ds}']

if pred_ds=='TRUB1_OE':
    base_mod = 'psi'
    other_mod = 'm6A'
    thresh_deltaS = 20
    df_deltaS = df_comp_pred[
        (df_comp_pred['name'] == base_mod)
        * (df_comp_pred['delta'] >= thresh_deltaS)
        ]
    ylim = 50
    xlabel = rf"Distance (bps) from ${{{dict_mod_display[base_mod]}}}$ site"\
             "\n"\
             rf"$\Delta S_{{{dict_mod_display[base_mod]}}}\geq{thresh_deltaS}$"
    fig_name = f'corr_func_deltaS_{other_mod}_from_{base_mod}_up_thresh{thresh_deltaS}.png'
elif pred_ds=='siTRUB1_input':
    base_mod = 'psi'
    other_mod = 'm6A'
    thresh_deltaS = 20
    df_deltaS = df_comp_pred[
        (df_comp_pred['name'] == base_mod)
        * (df_comp_pred['delta'] < -thresh_deltaS)
    ]
    ylim = 50
    xlabel = rf"Distance (bps) from ${{{dict_mod_display[base_mod]}}}$ site"\
             "\n"\
             rf"$\Delta S_{{{dict_mod_display[base_mod]}}}<-{thresh_deltaS}$"
    fig_name = f'corr_func_deltaS_{other_mod}_from_{base_mod}_down_thresh{thresh_deltaS}.png'
elif pred_ds=='METTL3_KO':
    base_mod = 'm6A'
    other_mod = 'psi'
    thresh_deltaS = 50
    df_deltaS = df_comp_pred[
        (df_comp_pred['name'] == base_mod)
        * (df_comp_pred['delta'] < -thresh_deltaS)
    ]
    ylim = 100
    xlabel = rf"Distance (bps) from ${{{dict_mod_display[base_mod]}}}$ site"\
             "\n"\
             rf"$\Delta S_{{{dict_mod_display[base_mod]}}}<-{thresh_deltaS}$"
    fig_name = f'corr_func_deltaS_{other_mod}_from_{base_mod}_down_thresh{thresh_deltaS}.png'
elif pred_ds=='METTL3_KD':
    base_mod = 'm6A'
    other_mod = 'psi'
    thresh_deltaS = 20
    df_deltaS = df_comp_pred[
        (df_comp_pred['name'] == base_mod)
        * (df_comp_pred['delta'] < -thresh_deltaS)
    ]
    ylim = 50
    xlabel = rf"Distance (bps) from ${{{dict_mod_display[base_mod]}}}$ site"\
             "\n"\
             rf"$\Delta S_{{{dict_mod_display[base_mod]}}}<-{thresh_deltaS}$"
    fig_name = f'corr_func_deltaS_{other_mod}_from_{base_mod}_down_thresh{thresh_deltaS}.png'



all_vec_r = []
all_vec_d = []
genes_deltaS = []
for _, this_row in tqdm(df_deltaS.iterrows()):
    this_chrom = this_row['chrom']
    this_chromStart = this_row['chromStart']
    this_chromEnd = this_row['chromEnd']
    this_strand = this_row['strand']

    sub_df_gene_bed = df_gene_bed[
        (df_gene_bed['chrom']==this_chrom)
        * (df_gene_bed['chromStart']<=this_chromStart)
        * (df_gene_bed['chromEnd']>=this_chromEnd)
        * (df_gene_bed['strand']==this_strand)
    ]

    if len(sub_df_gene_bed)==0:
        continue

    genes_deltaS.append(sub_df_gene_bed['name'].unique()[0])
    gene_chromStart = sub_df_gene_bed['chromStart'].min()
    gene_chromEnd = sub_df_gene_bed['chromEnd'].max()

    sub_df = df_comp_pred[
        (df_comp_pred['chrom']==this_chrom)
        * (df_comp_pred['name']==other_mod)
        * (df_comp_pred['strand']==this_strand)
        * (df_comp_pred['chromStart']>=gene_chromStart)
        * (df_comp_pred['chromStart']<gene_chromEnd)
    ]
    vec_r = sub_df['chromStart'].values - this_chromStart
    vec_d = sub_df['delta'].values
    if this_strand=='-':
        vec_r = -vec_r
    all_vec_r.append(vec_r)
    all_vec_d.append(vec_d)
all_vec_r = np.concatenate(all_vec_r)
all_vec_d = np.concatenate(all_vec_d)

max_dist = 4500
bin_width = 1000
bin_edges = np.arange(-max_dist, max_dist+bin_width, bin_width)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

binned_d = []
for this_bin_i in range(len(bin_edges)-1):
    bin_start = bin_edges[this_bin_i]
    bin_end = bin_edges[this_bin_i+1]
    binned_d.append(all_vec_d[(all_vec_r>=bin_start) * (all_vec_r<bin_end)])
binned_d_mean = [np.mean(this_bin) if len(this_bin) else np.nan for this_bin in binned_d]
binned_d_median = [np.median(this_bin) if len(this_bin) else np.nan for this_bin in binned_d]
binned_d = [this_bin if len(this_bin) else [0] for this_bin in binned_d]

plt.figure(figsize=(5, 5))
# plt.subplot(1, 2, 1)
plt.violinplot(binned_d, bin_centers, widths=int(bin_width/2))
plt.xlim([-max_dist, max_dist])
plt.ylim([-ylim, ylim])
plt.axhspan(-20, 20, color='gray', alpha=0.1)
# plt.axhline(y=25, c='gray')
# plt.axhline(y=-25, c='gray')
plt.xticks(bin_edges)
plt.xlabel(xlabel, fontsize=12)
plt.ylabel(f"$\Delta S_{{{dict_mod_display[other_mod]}}}$", fontsize=12)
plt.title(f"{pred_ds} cf. {comp_ds}", fontsize=15)
plt.savefig(os.path.join(img_out, fig_name), bbox_inches='tight')