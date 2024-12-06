import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

in_file = '/home/adrian/Data/UBS-Seq/41587_2023_2034_MOESM9_ESM.xlsx'
df_in = pd.read_excel(in_file, sheet_name='HEK293T-WT_sites', dtype={'chrom': str})
out_file = '/home/adrian/Data/UBS-Seq/HEK293T_WT.bed'

df_in = df_in[df_in['chrom'] != '.']

df_out = []
for _, row in df_in.iterrows():
    chrom = row['chrom']
    chromEnd = row['position']
    chromStart = chromEnd - 1
    strand = row['strand']
    modRatio = row['ratio'] * 100.0

    ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    if strand == '-':
        ref5mer = str(Seq(ref5mer).reverse_complement())
    # print(strand, ref5mer)

    df_out.append([chrom, chromStart, chromEnd, 'm5C', modRatio, strand, ref5mer])


df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False, float_format='%.1f')