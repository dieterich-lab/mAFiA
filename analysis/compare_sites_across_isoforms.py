import os
from glob import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)
from pybedtools import BedTool


thresh_delta = 20.0
thresh_confidence = 80.0
ds_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT'
gtf_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf'
gtf = BedTool(gtf_file)

def get_dict_gene_txId(transcript_dir):
    all_geneId_txId = [os.path.basename(path).rstrip('.bed').split('_') for path in glob(os.path.join(transcript_dir, '*.bed'))]
    unique_geneIds = list(set([pair[0] for pair in all_geneId_txId]))
    dict_gene_txId = {
        this_geneId: [pair[1] for pair in all_geneId_txId if pair[0]==this_geneId] for this_geneId in unique_geneIds
    }
    return dict_gene_txId


def load_dataframe(transcript_dir, geneId, txId, thresh_conf):
    df = pd.read_csv(os.path.join(transcript_dir, f'{geneId}_{txId}.bed'), sep='\t')
    df = df[df['confidence']>=thresh_confidence]
    # df.drop(columns='score', inplace=True)
    df.rename(columns={
        'coverage': f'coverage_{txId}',
        'modRatio': f'modRatio_{txId}',
        'confidence': f'confidence_{txId}'
    }, inplace=True)
    # df['transcript'] = txId
    # df['gene'] = geneId
    return df


def unpack_transcript_dataframe(in_dataframe):
    out_dataframe = []
    for _, this_row in in_dataframe.iterrows():
        unique_txIds = [k.split('_')[1] for k in this_row.keys() if k.split('_')[0] == 'modRatio']
        for this_txId in unique_txIds:
            if np.isnan(this_row[f'coverage_{this_txId}']):
                continue
            new_row = this_row[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer']]
            new_row['coverage'] = int(this_row[f'coverage_{this_txId}'])
            new_row['modRatio'] = this_row[f'modRatio_{this_txId}']
            new_row['confidence'] = this_row[f'confidence_{this_txId}']
            new_row['txId'] = this_txId
            new_row['geneId'] = this_row['geneId']
            out_dataframe.append(new_row)
    return pd.DataFrame(out_dataframe)


def get_diff_isoforms(transcript_dir, dict_gene_txId):
    diff_isoforms = []
    for this_geneId, this_geneTxIds in dict_gene_txId.items():
        if len(this_geneTxIds)>1:
            this_df = load_dataframe(transcript_dir, this_geneId, this_geneTxIds[0], thresh_confidence)
            for this_txId in this_geneTxIds[1:]:
                next_df = load_dataframe(transcript_dir, this_geneId, this_txId, thresh_confidence)
                this_df = pd.merge(this_df, next_df, how='outer',
                                   on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'])
                this_df[f'delta_{this_txId}_{this_geneTxIds[0]}'] = this_df[f'modRatio_{this_txId}'] - this_df[f'modRatio_{this_geneTxIds[0]}']

            abs_deltas = np.abs(this_df.filter(regex="delta_"))
            this_df['maxDelta'] = abs_deltas.max(axis=1)
            this_df_filtered = this_df[this_df['maxDelta']>=thresh_delta]

            if len(this_df_filtered):
                this_df_filtered['geneId'] = this_geneId
                diff_isoforms.append(unpack_transcript_dataframe(this_df_filtered))
            # print(this_df_filtered)
    return diff_isoforms

### todo ###
def annotate_site():
    return


chr_diff_isoforms = []
for this_chr in [str(c) for c in range(1, 23)] + ['X']:
    this_chr_transcript_dir = os.path.join(ds_dir, f'chr{this_chr}/transcript')
    this_chr_dict_gene_txId = get_dict_gene_txId(this_chr_transcript_dir)
    chr_diff_isoforms.extend(get_diff_isoforms(this_chr_transcript_dir, this_chr_dict_gene_txId))
chr_diff_isoforms = pd.concat(chr_diff_isoforms)

# for this_chr, this_chr_diff_isoforms in chr_diff_isoforms.items():
#     # print(this_chr)
#     for this_gene, this_gene_diff_isoforms in this_chr_diff_isoforms.items():
#         print('\n', this_gene)
#         print(this_gene_diff_isoforms)

