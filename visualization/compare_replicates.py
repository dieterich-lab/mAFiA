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

# chrom = '1'
chrom = 'ALL'
# train_ds = 'DRACH_v1'
# train_ds = 'ISA_retrain_GAACT_TGACT'

THRESH_CONF = 90

data_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/psico-mAFiA/HEK293'

df_wt = pd.read_csv(
    # f'{data_dir}/100_WT_0_IVT/chr{chrom}_confidence/mAFiA.sites.bed',
    f'{data_dir}/100_WT_0_IVT/chr{chrom}.mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str}
)
df_ivt = pd.read_csv(
    # f'{data_dir}/0_WT_100_IVT/chr{chrom}_confidence/mAFiA.sites.bed',
    f'{data_dir}/0_WT_100_IVT/chr{chrom}.mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str}
)
df_glori = pd.read_csv(
    f'/home/adrian/Data/GLORI/bed_files/GLORI.chr{chrom}.tsv',
    sep='\t',
    dtype={'chrom': str}
)
df_glori.rename(columns={'score': 'modRatio'}, inplace=True)

df_bid = pd.read_csv(
    f'/home/adrian/Data/BID_seq/BID_seq.bed',
    sep='\t',
    dtype={'chrom': str}
)
df_bid.rename(columns={'score': 'modRatio'}, inplace=True)

img_out = '/home/adrian/img_out/psico_mAFiA'
os.makedirs(img_out, exist_ok=True)

# six_motifs = [
#     'GGACT',
#     'GGACA',
#     'GAACT',
#     'AGACT',
#     'GGACC',
#     'TGACT'
# ]

# blacklist = [
#     'AGACC',
#     'TAACT',
#     'TAACA',
#     'AAACA',
#     'GAACC',
#     'AAACC'
# ]

blacklist = []
df_merged = pd.merge(df_wt, df_ivt, on=['chrom', 'chromStart', 'chromEnd', 'score', 'strand', 'ref5mer'], suffixes=['_wt', '_ivt'])
# df_merged_sel = df_merged[~df_merged['ref5mer'].isin(blacklist)]
df_merged_sel = df_merged[
    (df_merged['confidence_ivt']>=THRESH_CONF)
    * (df_merged['confidence_wt']>=THRESH_CONF)
    ]
#
# ### recalculate mod ratio by filtering out dubious pred5mers ###
#
# # flank_match_ratio = []
# corr_mod_ratio = []
# corr_coverage = []
# with pysam.AlignmentFile(bam_file, 'rb') as bam:
#     for _, row in tqdm(df_merged.iterrows()):
#         pred5mers = []
#         for this_col in bam.pileup(row['chrom'], row['chromStart'], row['chromEnd'], truncate=True):
#             if this_col.reference_pos==row['chromStart']:
#                 # total_reads = 0
#                 # corr_reads = 0
#                 col_mod_probs = []
#                 for pileupread in this_col.pileups:
#                     if not pileupread.is_del and not pileupread.is_refskip:
#                         # this_pred5mer = pileupread.alignment.get_forward_sequence()[pileupread.query_position-2:pileupread.query_position+3]
#                         this_pred5mer = pileupread.alignment.query_sequence[pileupread.query_position-2:pileupread.query_position+3]
#                         if pileupread.alignment.is_reverse:
#                             this_pred5mer = str(Seq(this_pred5mer).reverse_complement())
#
#                         if (
#                                 len(this_pred5mer)
#                                 and (pileupread.alignment.modified_bases is not None)
#                                 and len(pileupread.alignment.modified_bases)
#                                 and this_pred5mer==row['ref5mer']
#                                 # and (this_pred5mer[0]==row['ref5mer'][0])
#                                 # and (this_pred5mer[1]==row['ref5mer'][1])
#                                 # and (this_pred5mer[-2]==row['ref5mer'][-2])
#                                 # and (this_pred5mer[-1] == row['ref5mer'][-1])
#                         ):
#                             list_mod_probs = [pair[1] for pair in list(pileupread.alignment.modified_bases.values())[0] if pair[0]==pileupread.query_position]
#                             if len(list_mod_probs):
#                                 col_mod_probs.append(list_mod_probs[0])
#
#                         if len(this_pred5mer):
#                             pred5mers.append(this_pred5mer)
#
#                 print(row['ref5mer'], Counter(pred5mers))
#                 # print(col_mod_probs)
#                 if len(col_mod_probs):
#                     pos = sum([x>=204 for x in col_mod_probs])
#                     neg = sum([x<51 for x in col_mod_probs])
#                     corr_mod_ratio.append(np.round(pos / (pos+neg) * 100))
#                     corr_coverage.append(pos+neg)
#                     print(this_col.n, pos+neg)
#
#                     # corr_mod_ratio.append(np.round(np.mean([x>=128 for x in col_mod_probs])*100))
#                 else:
#                     corr_mod_ratio.append(0)
#                     corr_coverage.append(0)
#                         # total_bases += 1
#                         # if pileupread.alignment.get_forward_sequence()[pileupread.query_position]=='A':
#                         #     corr_bases += 1
#                 # basecall_accuracy.append(corr_bases/total_bases)
#                         # print('\tbase in read %s = %s' %
#                         #       (pileupread.alignment.query_name,
#                         #        pileupread.alignment.get_forward_sequence()[pileupread.query_position]))
#
#                 # num_all = len(pred5mers)
#                 # num_matched = len([x for x in pred5mers if
#                 #                    (x[0] == row['ref5mer'][0])
#                 #                    and (x[1] == row['ref5mer'][1])
#                 #                    and (x[-2] == row['ref5mer'][-2])
#                 #                    and (x[-1] == row['ref5mer'][-1])
#                 #                    ])
#                 # flank_match_ratio.append(num_matched/num_all)
#                 #
#                 # print(row['ref5mer'], this_col.n)
#                 # print(num_all, num_matched, f'{(num_matched/num_all):.2f}')
# df_merged['modRatio_ivt_corr'] = corr_mod_ratio
# df_merged['coverage_ivt_corr'] = corr_coverage

# plt.hist(flank_match_ratio, bins=50, range=[0, 1])

def scatter_plot(df_in, key_x, key_y, mod_type, fig_name):
    plt.figure(figsize=(5, 5))
    plt.plot(df_in[key_x], df_in[key_y], 'o', alpha=0.5)
    plt.xlim([-1, 101])
    plt.ylim([-1, 101])
    plt.xticks(np.linspace(0, 100, 5))
    plt.yticks(np.linspace(0, 100, 5))
    plt.xlabel("$S_{{{}}}$".format((key_x.split('_')[1]).upper()), fontsize=10)
    plt.ylabel("$S_{{{}}}$".format((key_y.split('_')[1]).upper()), fontsize=10)
    plt.suptitle(f'chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%', fontsize=12)
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
        plt.scatter(sub_df[key_x], sub_df[key_y], s=2, alpha=0.5)
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
    plt.suptitle(f'chr{chrom}, {len(df_in)} {mod_type} sites\nconf$\geq${THRESH_CONF}%', fontsize=12)
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

scatter_plot_by_motif(df_merged_sel, 'modRatio_wt', 'modRatio_ivt', 'm6A', m6A_motifs, 3, 6, f'm6A_WT_vs_IVT_conf{THRESH_CONF}.png')
scatter_plot_by_motif(df_merged_sel, 'modRatio_wt', 'modRatio_ivt', 'psi', psi_motifs, 2, 4, f'psi_WT_vs_IVT_conf{THRESH_CONF}.png')

### compare to GLORI ###
df_glori_wt = pd.merge(df_glori, df_wt, on=['chrom', 'chromStart', 'chromEnd', 'strand'], suffixes=['_glori', '_wt'])
df_glori_wt_sel = df_glori_wt[
    (df_glori_wt['confidence']>=THRESH_CONF)
]
scatter_plot_by_motif(df_glori_wt_sel, 'modRatio_glori', 'modRatio_wt', 'm6A', m6A_motifs, 3, 6, f'm6A_WT_vs_GLORI_conf{THRESH_CONF}.png')
scatter_plot(df_glori_wt_sel, 'modRatio_glori', 'modRatio_wt', 'm6A', f'm6A_WT_vs_GLORI_combined_conf{THRESH_CONF}.png')

### compare to Bid-Seq ###
df_bid_wt = pd.merge(df_bid, df_wt, on=['chrom', 'chromStart', 'chromEnd', 'strand', 'ref5mer'], suffixes=['_bid', '_wt'])
df_bid_wt_sel = df_bid_wt[
    (df_bid_wt['confidence']>=THRESH_CONF)
]
scatter_plot_by_motif(df_bid_wt_sel, 'modRatio_bid', 'modRatio_wt', 'psi', psi_motifs, 2, 4, f'psi_WT_vs_BID_conf{THRESH_CONF}.png')
scatter_plot(df_bid_wt_sel, 'modRatio_bid', 'modRatio_wt', 'psi', f'psi_WT_vs_BID_combined_conf{THRESH_CONF}.png')
