import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import pysam
import numpy as np
from tqdm import tqdm

df_wt = pd.read_csv(
    '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/100_WT_0_IVT/chr6/mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str})
df_ivt = pd.read_csv(
    '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/Mettl3-KO/chr6/mAFiA.sites.bed',
    sep='\t',
    dtype={'chrom': str})

bam_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH/100_WT_0_IVT/chr6/mAFiA.reads.bam'

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

basecall_accuracy = []
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for _, row in tqdm(df_merged.iterrows()):
        for this_col in bam.pileup(row['chrom'], row['chromStart'], row['chromEnd'], truncate=True):
            if this_col.reference_pos==row['chromStart']:
                total_bases = 0
                corr_bases = 0
                for pileupread in this_col.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        total_bases += 1
                        if pileupread.alignment.get_forward_sequence()[pileupread.query_position]=='A':
                            corr_bases += 1
                basecall_accuracy.append(corr_bases/total_bases)
                        # print('\tbase in read %s = %s' %
                        #       (pileupread.alignment.query_name,
                        #        pileupread.alignment.get_forward_sequence()[pileupread.query_position]))
df_merged['basecall_accuracy'] = basecall_accuracy

for motif in df_merged['ref5mer'].unique():
    df_merged_sel = df_merged[df_merged['ref5mer']==motif]
    df_merged_sel_acc = df_merged_sel[df_merged_sel['basecall_accuracy']<0.5]

    plt.figure(figsize=(5, 5))
    plt.scatter(df_merged_sel['modRatio_wt'], df_merged_sel['modRatio_ivt'], alpha=0.5)
    plt.scatter(df_merged_sel_acc['modRatio_wt'], df_merged_sel_acc['modRatio_ivt'], c='r', alpha=0.5)
    plt.xlim([-1, 101])
    plt.ylim([-1, 101])
    plt.xlabel('$S_{WT}$')
    plt.ylabel('$S_{IVT}$')
    plt.title(motif)
    # plt.title('chr1 DRACH sites')

print(Counter(df_merged[(df_merged['modRatio_ivt']>50)*(df_merged['modRatio_wt']<50)]['ref5mer']))

# sub_df = df_merged_sel[df_merged_sel['basecall_accuracy']>=0.5]
# plt.figure(figsize=(5, 5))
# plt.scatter(sub_df['modRatio_wt'].values, sub_df['modRatio_ivt'].values)
# plt.xlim([-1, 101])
# plt.ylim([-1, 101])
