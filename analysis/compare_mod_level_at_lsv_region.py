import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

thresh_conf = 60.0

day = 'day21'
majiq_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/majiq/TAC/voila_modulize_{day}/summary.tsv'
sham_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_{day}/chrALL.mAFiA.sites.bed'
tac_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_{day}/chrALL.mAFiA.sites.bed'

df_majiq = pd.read_csv(majiq_file, sep='\t', comment='#')
df_sham = pd.read_csv(sham_file, sep='\t', dtype={'chrom': str})
df_tac = pd.read_csv(tac_file, sep='\t', dtype={'chrom': str})

def filter_df(in_df, chrom, strand, chromStart, chromEnd):
    return in_df[
        (in_df['chrom'] == chrom)
        * (in_df['strand'] == strand)
        * (in_df['chromStart'] >= chromStart)
        * (in_df['chromEnd'] <= chromEnd)
        * (in_df['confidence'] >= thresh_conf)
    ]

for _, row in df_majiq.iterrows():
    chrom = row['seqid']
    strand = row['strand']
    lsv_regions = row['lsv_id'].split(';')
    for this_region in lsv_regions:
        source_target = this_region.split(':')[-2]
        chromStart, chromEnd = [int(x) for x in this_region.split(':')[-1].split('-')]
        chromStart -= 1
        chromEnd -= 1

        sub_df_sham = filter_df(df_sham, chrom, strand, chromStart, chromEnd)
        sub_df_tac = filter_df(df_tac, chrom, strand, chromStart, chromEnd)
        merged_df_sham_tac = pd.merge(sub_df_sham, sub_df_tac,
                                      on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                                      suffixes=['_sham', '_tac']
                                      )
        merged_df_sham_tac['delta'] = merged_df_sham_tac['modRatio_tac'] - merged_df_sham_tac['modRatio_sham']

        if len(merged_df_sham_tac)==0:
            continue

        plt.figure(figsize=(10, 5))
        for this_subplot, this_mod in enumerate(['m6A', 'psi']):
            plt.subplot(1, 2, this_subplot+1)
            this_df = merged_df_sham_tac[merged_df_sham_tac['name']==this_mod]
            plt.plot(this_df['delta'], this_df['modRatio_sham'], '.')
            plt.xlim([-20, 20])
            plt.ylim([-5, 100])
            plt.xlabel('${\Delta}S$')
            plt.ylabel('$S_{SHAM}$')
            plt.title(this_mod)
        plt.suptitle(f"{this_region}\n{row['module_event_combination']}", fontsize=15)
