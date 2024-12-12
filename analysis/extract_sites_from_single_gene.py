import os
import pandas as pd
import re
from functools import reduce
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'frame',
    'description'
]
df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

mafia_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

sel_gene = 'S100a1'
df_gene_sel = df_gene[df_gene['gene'] == sel_gene]
sel_chrom, sel_chromStart, sel_chromEnd, sel_strand = df_gene_sel[['chrom', 'chromStart', 'chromEnd', 'strand']].values[0]

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
conditions = {
    'TAC': ['SHAM', 'TAC'],
    'HFpEF': ['ctrl', 'HFpEF']
}
condition_suffix = {
    'TAC': {
        'SHAM': 'HFrEF_ctrl',
        'TAC': 'HFrEF_disease'
    },
    'HFpEF': {
        'ctrl': 'HFpEF_ctrl',
        'HFpEF': 'HFpEF_disease'
    }
}

all_dfs = []
for this_ds in conditions.keys():
    for this_cond in conditions[this_ds]:
        this_filepath = os.path.join(res_dir, this_ds, f'{this_cond}_merged', 'chrALL.mAFiA.sites.bed')
        this_df = pd.read_csv(this_filepath, sep='\t', dtype={'chrom': str})
        this_df_sel_gene = this_df[
            (this_df['chrom'] == sel_chrom)
            * (this_df['chromStart'] >= sel_chromStart)
            * (this_df['chromEnd'] <= sel_chromEnd)
            * (this_df['strand'] == sel_strand)
        ]
        this_df_sel_gene.rename(columns={
            'coverage': f'coverage_{condition_suffix[this_ds][this_cond]}',
            'modRatio': f'modRatio_{condition_suffix[this_ds][this_cond]}',
            'confidence': f'confidence_{condition_suffix[this_ds][this_cond]}',
        }, inplace=True)
        all_dfs.append(this_df_sel_gene)
df_merged = df_merged = reduce(lambda left, right: pd.merge(left, right, on=mafia_fields, how='outer'), all_dfs)

df_merged.to_csv(os.path.join(res_dir, f'mod_sites_{sel_gene}.tsv'), sep='\t', index=False)


