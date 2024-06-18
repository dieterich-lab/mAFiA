import os
import pandas as pd
pd.set_option('display.max_columns', None)
from pybedtools import BedTool
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np


dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

# gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
# gtf = BedTool(gtf_file)


thresh_conf = 50.0
day = 'day56'

majiq_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/majiq/TAC/voila_modulize_{day}/summary.tsv'
sham_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_{day}/chrALL.mAFiA.sites.bed'
tac_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_{day}/chrALL.mAFiA.sites.bed'
img_out = f'/home/adrian/img_out/TAC_splice_site_mod_changes/{day}'
os.makedirs(img_out, exist_ok=True)

df_majiq = pd.read_csv(majiq_file, sep='\t', comment='#')
df_sham = pd.read_csv(sham_file, sep='\t', dtype={'chrom': str})
df_tac = pd.read_csv(tac_file, sep='\t', dtype={'chrom': str})

def filter_df(in_df, chrom, strand, chromStart, chromEnd, thresh_conf):
    return in_df[
        (in_df['chrom'] == chrom)
        * (in_df['strand'] == strand)
        * (in_df['chromStart'] >= chromStart)
        * (in_df['chromEnd'] <= chromEnd)
        * (in_df['confidence'] >= thresh_conf)
    ]


for _, row in df_majiq.iterrows():
    gene_name = row['gene_name']
    chrom = row['seqid']
    strand = row['strand']
    lsv_regions = row['lsv_id'].split(';')
    for this_lsv_region in lsv_regions:
        source_target = this_lsv_region.split(':')[-2]
        chromStart, chromEnd = [int(x) for x in this_lsv_region.split(':')[-1].split('-')]
        chromStart -= 1
        chromEnd -= 1

        sub_df_sham = filter_df(df_sham, chrom, strand, chromStart, chromEnd, thresh_conf)
        sub_df_tac = filter_df(df_tac, chrom, strand, chromStart, chromEnd, thresh_conf)
        merged_df_sham_tac = pd.merge(sub_df_sham, sub_df_tac,
                                      on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                                      how='outer',
                                      suffixes=['_sham', '_tac']
                                      )
        merged_df_sham_tac['delta'] = merged_df_sham_tac['modRatio_tac'] - merged_df_sham_tac['modRatio_sham']

        if len(merged_df_sham_tac)==0:
            continue

        # plt.figure(figsize=(10, 5))
        # for this_subplot, this_mod in enumerate(['m6A', 'psi']):
        #     plt.subplot(1, 2, this_subplot+1)
        #     this_df = merged_df_sham_tac[merged_df_sham_tac['name']==this_mod]
        #     plt.plot(this_df['delta'], this_df['modRatio_sham'], '.')
        #     plt.xlim([-20, 20])
        #     plt.ylim([-5, 100])
        #     plt.xlabel('${\Delta}S$')
        #     plt.ylabel('$S_{SHAM}$')
        #     plt.title(this_mod)
        # plt.suptitle(f"{this_lsv_region}\n{row['module_event_combination']}", fontsize=15)

        # df_bed = pd.DataFrame.from_dict({
        #     'chrom': [chrom],
        #     'chromStart': [chromStart],
        #     'chromEnd': [chromEnd],
        #     'name': [this_lsv_region],
        #     'score': ['.'],
        #     'strand': [strand]
        # })
        # annotations = gtf.intersect(BedTool.from_dataframe(df_bed))
        # exons = [annot for annot in annotations if
        #         ('exon_number' in annot.attrs.keys()) and ('exon_id' in annot.attrs.keys())]
        # unique_tx_ids = list(set([this_exon.attrs['transcript_id'] for this_exon in exons]))
        # tx_colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_tx_ids)))

        plt.figure(figsize=(10, 10))
        for this_subplot, this_mod in enumerate(['m6A', 'psi']):
            plt.subplot(2, 1, this_subplot+1)
            this_df = merged_df_sham_tac[merged_df_sham_tac['name']==this_mod]
            labelled = False
            for _, this_site in this_df.iterrows():
                x = this_site['chromEnd']
                if not labelled:
                    plt.plot(x, this_site['modRatio_sham'], 'b+', label=f'SHAM_{day}')
                    plt.plot(x, this_site['modRatio_tac'], 'r+', label=f'TAC_{day}')
                    labelled = True
                else:
                    plt.plot(x, this_site['modRatio_sham'], 'b+')
                    plt.plot(x, this_site['modRatio_tac'], 'r+')
                lc = 'r' if this_site['delta']>=0 else 'b'
                plt.plot((x, x), (this_site['modRatio_sham'], this_site['modRatio_tac']), '-', c=lc)
            # plt.xlim([-20, 20])
            # used_exon_ids = []

            # for i, this_tx_id in enumerate(unique_tx_ids):
            #     tx_exons = [this_exon for this_exon in exons if this_exon.attrs['transcript_id']==this_tx_id]
            #     for this_exon in tx_exons:
            #     # if this_exon.attrs['exon_id'] not in used_exon_ids:
            #         plt.axvline(this_exon.start, c=tx_colors[i])
            #         plt.axvline(this_exon.end, c=tx_colors[i])
            #         # plt.text(0.5*(this_exon.start+this_exon.end), 100, this_exon.attrs['exon_number'])
            #     # used_exon_ids.append(this_exon.attrs['exon_id'])

            plt.ylim([-5, 100])
            plt.xlabel('Exon position', fontsize=12)
            plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
            plt.axvline(chromStart, c='gray')
            plt.axvline(chromEnd, c='gray')
            if strand == '+':
                plt.xticks([chromStart, chromEnd], [f"5'\nchr{chrom}: {chromStart+1}", f"3'\nchr{chrom}: {chromEnd+1}"])
                plt.legend(loc='upper left')
            else:
                plt.xticks([chromStart, chromEnd], [f"3'\nchr{chrom}: {chromStart+1}", f"5'\nchr{chrom}: {chromEnd+1}"])
                plt.legend(loc='upper right')
            region_len = chromEnd - chromStart
            plt.xlim([chromStart-0.1*region_len, chromEnd+0.1*region_len])

        plt.suptitle(f"{gene_name}\n{this_lsv_region}\n{row['module_event_combination']}", fontsize=15)
        plt.savefig(os.path.join(img_out, f'{gene_name}_{this_lsv_region}.png'), bbox_inches='tight')
        plt.close('all')
