import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import pysam
import numpy as np
from tqdm import tqdm
import os
from Bio.Seq import Seq

chrom = '6'
train_ds = 'DRACH_v1'
# train_ds = 'ISA_retrain_GAACT_TGACT'

df_wt = pd.read_csv(
    # '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/100_WT_0_IVT/chr1/mAFiA.sites.bed',
    f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/100_WT_0_IVT/chr{chrom}/mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str})
df_ivt = pd.read_csv(
    # '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/GLORI/0_WT_100_IVT/chr1/mAFiA.sites.bed',
    # f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/{train_ds}/0_WT_100_IVT/chr{chrom}/mAFiA.sites.bed',
    f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/{train_ds}/0_WT_100_IVT/chr{chrom}/mAFiA.sites.bed',
    # '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/Mettl3-KO/chr6/mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str})

bam_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/{train_ds}/0_WT_100_IVT/chr{chrom}/mAFiA.reads.bam'

img_out = '/home/adrian/img_out/compare_IVT'
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
# df_merged_sel = df_merged

### recalculate mod ratio by filtering out dubious pred5mers ###

# flank_match_ratio = []
corr_mod_ratio = []
corr_coverage = []
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for _, row in tqdm(df_merged.iterrows()):
        pred5mers = []
        for this_col in bam.pileup(row['chrom'], row['chromStart'], row['chromEnd'], truncate=True):
            if this_col.reference_pos==row['chromStart']:
                # total_reads = 0
                # corr_reads = 0
                col_mod_probs = []
                for pileupread in this_col.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # this_pred5mer = pileupread.alignment.get_forward_sequence()[pileupread.query_position-2:pileupread.query_position+3]
                        this_pred5mer = pileupread.alignment.query_sequence[pileupread.query_position-2:pileupread.query_position+3]
                        if pileupread.alignment.is_reverse:
                            this_pred5mer = str(Seq(this_pred5mer).reverse_complement())

                        if (
                                len(this_pred5mer)
                                and (pileupread.alignment.modified_bases is not None)
                                and len(pileupread.alignment.modified_bases)
                                and this_pred5mer==row['ref5mer']
                                # and (this_pred5mer[0]==row['ref5mer'][0])
                                # and (this_pred5mer[1]==row['ref5mer'][1])
                                # and (this_pred5mer[-2]==row['ref5mer'][-2])
                                # and (this_pred5mer[-1] == row['ref5mer'][-1])
                        ):
                            list_mod_probs = [pair[1] for pair in list(pileupread.alignment.modified_bases.values())[0] if pair[0]==pileupread.query_position]
                            if len(list_mod_probs):
                                col_mod_probs.append(list_mod_probs[0])

                        if len(this_pred5mer):
                            pred5mers.append(this_pred5mer)

                print(row['ref5mer'], Counter(pred5mers))
                # print(col_mod_probs)
                if len(col_mod_probs):
                    pos = sum([x>=204 for x in col_mod_probs])
                    neg = sum([x<51 for x in col_mod_probs])
                    corr_mod_ratio.append(np.round(pos / (pos+neg) * 100))
                    corr_coverage.append(pos+neg)
                    print(this_col.n, pos+neg)

                    # corr_mod_ratio.append(np.round(np.mean([x>=128 for x in col_mod_probs])*100))
                else:
                    corr_mod_ratio.append(0)
                    corr_coverage.append(0)
                        # total_bases += 1
                        # if pileupread.alignment.get_forward_sequence()[pileupread.query_position]=='A':
                        #     corr_bases += 1
                # basecall_accuracy.append(corr_bases/total_bases)
                        # print('\tbase in read %s = %s' %
                        #       (pileupread.alignment.query_name,
                        #        pileupread.alignment.get_forward_sequence()[pileupread.query_position]))

                # num_all = len(pred5mers)
                # num_matched = len([x for x in pred5mers if
                #                    (x[0] == row['ref5mer'][0])
                #                    and (x[1] == row['ref5mer'][1])
                #                    and (x[-2] == row['ref5mer'][-2])
                #                    and (x[-1] == row['ref5mer'][-1])
                #                    ])
                # flank_match_ratio.append(num_matched/num_all)
                #
                # print(row['ref5mer'], this_col.n)
                # print(num_all, num_matched, f'{(num_matched/num_all):.2f}')
df_merged['modRatio_ivt_corr'] = corr_mod_ratio
df_merged['coverage_ivt_corr'] = corr_coverage

# plt.hist(flank_match_ratio, bins=50, range=[0, 1])

ordered_motifs = [
    'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
    'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
     'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
]
# ordered_motifs = ordered_motifs + [pair[0] for pair in Counter(df_wt['ref5mer']).most_common() if pair[0] not in ordered_motifs]

num_row = 3
num_col = 6

thresh_err = 25
motif_counts = Counter(df_merged['ref5mer'])
motif_err_rate = {motif: (df_merged[df_merged['ref5mer']==motif]['modRatio_ivt']>=thresh_err).sum() / motif_counts[motif] for motif in motif_counts.keys()}

plt.figure(figsize=(12, 6))
for ind, motif in enumerate(ordered_motifs):
    df_merged_sel = df_merged[df_merged['ref5mer']==motif]
    df_merged_sel_acc = df_merged_sel[df_merged_sel['coverage_ivt_corr']>=50]

    plt.subplot(num_row, num_col, ind+1)
    # plt.scatter(df_merged_sel['modRatio_wt'], df_merged_sel['modRatio_ivt_corr'], s=1.5, alpha=0.5)
    plt.scatter(df_merged_sel_acc['modRatio_wt'], df_merged_sel_acc['modRatio_ivt'], s=0.5, c='r', alpha=0.5)
    plt.text(1, 90, f"{motif.replace('T', 'U')} ({motif_err_rate[motif]:.2f})", fontsize=10)
    plt.axhline(y=thresh_err, linestyle='--', c='r')
    # plt.legend(loc='upper left', handlelength=0.1, fontsize=10)
    plt.xlim([-1, 101])
    plt.ylim([-1, 101])
    plt.xticks(np.linspace(0, 100, 5))
    plt.yticks(np.linspace(0, 100, 5))
    if ind >= ((num_row-1)*num_col):
        plt.xlabel('$S_{WT}$', fontsize=10)
    if ind % num_col==0:
        plt.ylabel('$S_{IVT}$', fontsize=10)
    # plt.title(motif)
    # plt.title('chr1 DRACH sites')
plt.suptitle(f'chr{chrom} DRACH sites', fontsize=12)
plt.savefig(os.path.join(img_out, f'scatter_WT_vs_IVT_{train_ds}.png'), bbox_inches='tight')

# print(Counter(df_merged[(df_merged['modRatio_ivt']>50)]['ref5mer']))
# print(Counter(df_merged[(df_merged['modRatio_ivt']>50)*(df_merged['modRatio_wt']<50)]['ref5mer']))

# sub_df = df_merged_sel[df_merged_sel['basecall_accuracy']>=0.5]
# plt.figure(figsize=(5, 5))
# plt.scatter(sub_df['modRatio_wt'].values, sub_df['modRatio_ivt'].values)
# plt.xlim([-1, 101])
# plt.ylim([-1, 101])
