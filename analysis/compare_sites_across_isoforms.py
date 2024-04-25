import os
from glob import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)
from pybedtools import BedTool
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


thresh_delta = 10.0
thresh_confidence = 50.0
ds_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT'
gtf_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf'
gtf = BedTool(gtf_file)

annotated_df_file = os.path.join(ds_dir, 'diff_isoforms_annotated.tsv')

img_out = '/home/adrian/img_out/transcript_mod_profile'
os.makedirs(img_out, exist_ok=True)

def get_dict_gene_txId(transcript_dir):
    all_geneId_txId = [os.path.basename(path).rstrip('.bed').split('_') for path in glob(os.path.join(transcript_dir, '*.bed'))]
    unique_geneIds = list(set([pair[0] for pair in all_geneId_txId]))
    dict_gene_txId = {
        this_geneId: [pair[1] for pair in all_geneId_txId if pair[0]==this_geneId] for this_geneId in unique_geneIds
    }
    return dict_gene_txId


def load_dataframe(transcript_dir, geneId, txId, thresh_conf):
    df = pd.read_csv(os.path.join(transcript_dir, f'{geneId}_{txId}.bed'), sep='\t', dtype={'chrom': str})
    df = df[df['confidence']>=thresh_conf]
    # df.drop(columns='score', inplace=True)
    df.rename(columns={
        'coverage': f'coverage_{txId}',
        'modRatio': f'modRatio_{txId}',
        'confidence': f'confidence_{txId}'
    }, inplace=True)
    # df['transcript'] = txId
    # df['gene'] = geneId
    return df


def unpack_transcript_dataframe(in_dataframe, geneId):
    out_dataframe = []
    for _, this_row in in_dataframe.iterrows():
        unique_txIds = [k.split('_')[1] for k in this_row.keys() if k.split('_')[0] == 'modRatio']
        for this_txId in unique_txIds:
            if np.isnan(this_row[f'coverage_{this_txId}']):
                continue
            new_row = this_row[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer']].copy()
            new_row['coverage'] = int(this_row[f'coverage_{this_txId}'])
            new_row['modRatio'] = this_row[f'modRatio_{this_txId}']
            new_row['confidence'] = this_row[f'confidence_{this_txId}']
            new_row['txId'] = this_txId
            new_row['geneId'] = geneId
            out_dataframe.append(new_row)
    return pd.DataFrame(out_dataframe)


def get_diff_isoforms(transcript_dir, dict_gene_txId):
    print('Collecting differential isoforms...')
    diff_isoforms = []
    for this_geneId, this_geneTxIds in tqdm(dict_gene_txId.items()):
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
                # this_df_filtered['geneId'] = this_geneId
                diff_isoforms.append(unpack_transcript_dataframe(this_df_filtered, this_geneId))
    if len(diff_isoforms):
        diff_isoforms = pd.concat(diff_isoforms).sort_values(['chromStart', 'txId'])
    return diff_isoforms


def annotate_sites(in_df):
    print('Annotating sites...')

    gtf = BedTool(gtf_file)
    unique_chromStarts = in_df['chromStart'].unique()

    out_df = []
    for this_chromStart in tqdm(unique_chromStarts):
        sub_df = in_df[in_df['chromStart']==this_chromStart]
        this_site_bed = BedTool.from_dataframe(sub_df)
        annotations = gtf.intersect(this_site_bed)
        for _, this_row in sub_df.iterrows():
            this_txId = this_row['txId']
            this_tx_annots = [annot for annot in annotations if annot.attrs.get('transcript_id')==this_txId]
            this_location = [annot.fields[2] for annot in this_tx_annots if annot.fields[2] not in ['transcript', 'exon']]
            this_biotype = [annot.attrs.get('transcript_biotype') for annot in this_tx_annots]
            this_row['location'] = ', '.join(list(set(this_location)))
            this_row['biotype'] = ', '.join(list(set(this_biotype)))
            out_df.append(this_row)
    out_df = pd.DataFrame(out_df)
    return out_df


def get_exon_stoichiometry(pos, modRatio, ranges):
    exon_pos = []
    exon_modRatio = []
    for this_range in ranges:
        mask = (pos>=this_range[0]) * (pos<=this_range[1])
        exon_pos.append(np.mean(this_range))
        if np.sum(mask)==0:
            exon_modRatio.append(0.0)
        else:
            # exon_pos.append(np.mean(pos[mask]))
            exon_modRatio.append(np.mean(modRatio[mask]))
    return exon_pos, exon_modRatio


def plot_transcripts(transcript_dir, geneId, list_txIds, strand):
    mod_colors = {
        'm6A': 'r',
        'psi': 'g'
    }

    # transcript_dir = os.path.join(ds_dir, 'chr1/transcript')
    # geneId = 'ENSG00000143774'
    # list_txIds = ['ENST00000312726', 'ENST00000391865', 'ENST00000471270']

    gene_annotations = [annot for annot in gtf.filter(lambda x: x.attrs['gene_id']==geneId)]
    gene_interval = [annot for annot in gene_annotations if annot.fields[2]=='gene'][0]
    gene_name = gene_interval.attrs['gene_name']
    if strand=='+':
        gene_start = gene_interval.start
        gene_end = gene_interval.end
    else:
        gene_start = gene_interval.end
        gene_end = gene_interval.start

    num_rows = len(list_txIds)

    fig = plt.figure(figsize=(10, num_rows*4))
    for ind, this_txId in enumerate(list_txIds):
        tx_annotations = [annot for annot in gene_annotations if annot.attrs.get('transcript_id') == this_txId]
        tx_biotype = [annot.attrs['transcript_biotype'] for annot in tx_annotations if annot.fields[2]=='transcript'][0]
        exon_ranges = [(annot.start, annot.end) for annot in tx_annotations if annot.fields[2]=='exon']
        cds_ranges = [(annot.start, annot.end) for annot in tx_annotations if annot.fields[2]=='CDS']
        utr_ranges = [(annot.start, annot.end) for annot in tx_annotations if annot.fields[2] in ['five_prime_utr', 'three_prime_utr']]

        if strand=='-':
            exon_ranges = [this_range[::-1] for this_range in exon_ranges]
            cds_ranges = [this_range[::-1] for this_range in cds_ranges]
            utr_ranges = [this_range[::-1] for this_range in utr_ranges]

        plt.subplot(num_rows, 1, ind+1)
        this_tx_df = load_dataframe(transcript_dir, geneId, this_txId, thresh_conf=thresh_confidence)
        for mod in ['m6A', 'psi']:
            vec_pos, vec_modRatio = this_tx_df[this_tx_df['name']==mod][['chromEnd', f'modRatio_{this_txId}']].values.T
            if strand=='-':
                vec_pos = vec_pos[::-1]
                vec_modRatio = vec_modRatio[::-1]
            # exon_pos, exon_modRatio = get_exon_stoichiometry(vec_pos, vec_modRatio, exon_ranges)
            plt.plot(vec_pos, vec_modRatio, f'{mod_colors[mod]}.-', label=mod)
            for exon in exon_ranges:
                plt.axvspan(xmin=exon[0], xmax=exon[1], color='grey', alpha=0.1)
            for cds in cds_ranges:
                plt.axvspan(xmin=cds[0], xmax=cds[1], color='blue', alpha=0.1)
            for utr in utr_ranges:
                plt.axvspan(xmin=utr[0], xmax=utr[1], color='yellow', alpha=0.1)

        plt.title(f'{this_txId} ({tx_biotype})', fontsize=15)
        plt.xlim([gene_start, gene_end])
        plt.ylim([-5, 105])
        if ind==(num_rows-1):
            plt.xlabel('Genomic coord.', fontsize=12)
        else:
            plt.xticks([])
        plt.ylabel('Mod. Ratio', fontsize=12)
        plt.legend(loc='upper left')
    plt.suptitle(f'{geneId}: {gene_name}', fontsize=15)

    return fig


if os.path.exists(annotated_df_file):
    chr_diff_isoforms = pd.read_csv(annotated_df_file, sep='\t', dtype={'chrom': str})
else:
    chr_diff_isoforms = []
    for this_chr in [str(c) for c in range(1, 23)] + ['X']:
        print(f'\nchr{this_chr}')
        this_chr_transcript_dir = os.path.join(ds_dir, f'chr{this_chr}/transcript')
        this_chr_dict_gene_txId = get_dict_gene_txId(this_chr_transcript_dir)
        diff_isoforms = get_diff_isoforms(this_chr_transcript_dir, this_chr_dict_gene_txId)
        if len(diff_isoforms):
            annotated_isoforms = annotate_sites(diff_isoforms)
            chr_diff_isoforms.append(annotated_isoforms)
    chr_diff_isoforms = pd.concat(chr_diff_isoforms)
    chr_diff_isoforms.to_csv(annotated_df_file, sep='\t', index=False)

sel_geneIds = list(chr_diff_isoforms['geneId'].unique())
# sel_geneIds = list(chr_diff_isoforms[chr_diff_isoforms['biotype']=='retained_intron']['geneId'].unique())
sel_geneIds.sort()

for this_geneId in sel_geneIds:
    sub_df = chr_diff_isoforms[chr_diff_isoforms['geneId']==this_geneId]
    this_list_txIds = list(sub_df['txId'].unique())
    this_chr = sub_df['chrom'].values[0]
    this_transcript_dir = os.path.join(ds_dir, f'chr{this_chr}/transcript')
    this_strand = sub_df['strand'].values[0]
    this_fig = plot_transcripts(this_transcript_dir, this_geneId, this_list_txIds, this_strand)
    this_fig.savefig(os.path.join(img_out, f'{this_geneId}.png'), bbox_inches='tight')
    plt.close('all')