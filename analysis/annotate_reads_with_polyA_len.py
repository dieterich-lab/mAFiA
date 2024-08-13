import pandas as pd
import pysam
from tqdm import tqdm
import re
import os

ds = 'TAC'
cond = 'SHAM_day56'
polyA_dir = f"/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/polyA"
polyA_file = os.path.join(polyA_dir, f"{cond}_read2polyA_length.txt")
bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/{cond}/chrALL.mAFiA.reads.bam'
gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'

out_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA'
os.makedirs(out_dir, exist_ok=True)

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

df_polyA = pd.read_csv(polyA_file, sep='\t', names=['read_id', 'polyA_length'])
read_ids_with_polyA = list(df_polyA['read_id'].values)

df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

gene_to_read = {}
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for _, this_row in tqdm(df_gene.iterrows()):
        if this_row['strand'] == '+':
            flag_required = 0
        else:
            flag_required = 16
        collected_read_ids = []
        for this_read in bam.fetch(this_row['chrom'], this_row['chromStart'], this_row['chromEnd']):
            if this_read.flag == flag_required:
                collected_read_ids.append(this_read.query_name)
        if len(collected_read_ids):
            gene_to_read[this_row['gene']] = collected_read_ids

read_to_gene = {
    read: gene
    for gene, reads in gene_to_read.items()
    for read in reads
}

df_polyA['gene'] = [read_to_gene.get(this_read) for this_read in df_polyA['read_id']]
df_polyA_sel = df_polyA[~df_polyA['gene'].isna()]

df_polyA_sel.to_csv(os.path.join(out_dir, f'polyA_reads_annotated_{ds}_{cond}.tsv'), sep='\t', index=False)
print(f'{len(df_polyA_sel)} reads written')
