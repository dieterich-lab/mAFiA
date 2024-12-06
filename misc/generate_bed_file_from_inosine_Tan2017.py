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

in_file = '/home/adrian/Data/inosine/Tan2017_Supplementary_File_7_GRCh38.tsv'
df_in = pd.read_csv(in_file, sep='\t', dtype={'chrom': str})
out_file = '/home/adrian/Data/inosine/HEK293_Tan2017_OE.bed'

df_in['chrom'] = [this_chrom.lstrip('chr') for this_chrom in df_in['chrom']]

df_out = []
for _, row in df_in.iterrows():
    chrom = row['chrom']
    chromStart = row['chromStart']
    chromEnd = row['chromEnd']
    modRatio = row['OE_EditingLevel'] * 100.0

    ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    if ref5mer[2] == 'A':
        strand = '+'
    elif ref5mer[2] == 'T':
        strand = '-'
        ref5mer = str(Seq(ref5mer).reverse_complement())

    # print(strand, ref5mer)

    df_out.append([chrom, chromStart, chromEnd, 'ino', modRatio, strand, ref5mer])


df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False, float_format='%.1f')