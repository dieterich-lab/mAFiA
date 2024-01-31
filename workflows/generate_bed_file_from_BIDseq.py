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

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

in_file = '/home/adrian/Data/BID_seq/GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx'
df_in = pd.read_excel(in_file, skiprows=range(3), usecols=['chr', 'pos', 'strand', 'Motif_1', 'Frac_Ave %'])

out_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/psU/site_annotations/BID_seq.bed'

num_false = 0
df_out = []
for _, row in df_in.iterrows():
    chrom = row['chr'].lstrip('chr')
    chromEnd = row['pos']
    chromStart = chromEnd - 1
    strand = row['strand']

    # ref5mer = ref[chrom][chromStart-2:chromStart+3]
    # if strand=='-':
    #     ref5mer = str(Seq(ref5mer).reverse_complement())

    ref9mer = ref[chrom][chromStart-4:chromStart+5]
    if strand=='-':
        ref9mer = str(Seq(ref9mer).reverse_complement())

    span = re.search(row['Motif_1'], ref9mer).span()
    if span[0]!=2:
        continue
        # shift = span[0] - 2
        # if strand=='+':
        #     chromStart += shift
        #     chromEnd += shift
        # else:
        #     chromStart -= shift
        #     chromEnd -= shift
    ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    if strand == '-':
        ref5mer = str(Seq(ref5mer).reverse_complement())

    # print(ref5mer==row['Motif_1'])

    df_out.append([chrom, chromStart, chromEnd, 'psU', round(row['Frac_Ave %'], 1), strand, ref5mer])
df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False)