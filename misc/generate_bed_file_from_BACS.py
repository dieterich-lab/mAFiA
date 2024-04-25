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

in_file = '/home/adrian/Data/BACS/mRNA_psi_list.xlsx'
out_file = '/home/adrian/Data/BACS/BACS_HeLa_WT.bed'
ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'

ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

df_in = pd.read_excel(in_file, skiprows=range(2), usecols=['chr', 'pos', 'strand', 'motif', 'Ψ level'])

df_out = []
for _, row in df_in.iterrows():
    chrom = row['chr'].lstrip('chr')
    chromEnd = row['pos']
    chromStart = chromEnd - 1
    strand = row['strand']

    ref5mer = ref[chrom][chromStart-2:chromStart+3]
    if strand=='-':
        ref5mer = str(Seq(ref5mer).reverse_complement())

    if ref5mer==row['motif']:
        # print(ref5mer, row['motif'])
        df_out.append([chrom, chromStart, chromEnd, 'psi', round(row['Ψ level']*100.0, 1), strand, ref5mer])
df_out = pd.DataFrame(df_out, columns=bed_fields)

df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False)