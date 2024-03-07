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

# ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref_file = '/home/adrian/Data/GRCm38_102/GRCm38_102.fa'

ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

# in_file = '/home/adrian/Data/BID_seq/GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx'
in_file = '/home/adrian/Data/BID_seq/41587_2022_1505_MOESM3_ESM.xlsx'
sheet_name = 'Supplementary Table 15'
out_file = '/home/adrian/Data/BID_seq/BID_seq_mouse_heart.bed'

df_in = pd.read_excel(in_file, sheet_name=sheet_name, skiprows=range(3), usecols=['chr', 'pos', 'strand', 'Motif_1', 'Motif_2', 'Frac_Ave %'])
df_in = df_in[df_in['Motif_2'].isna()]

def expanded_search(this_row, this_ref, test_motif):
    chrom = this_row['chr'].lstrip('chr')
    chromEnd = this_row['pos']
    chromStart = chromEnd - 1
    strand = this_row['strand']
    ref7mer = this_ref[chrom][chromStart-3:chromStart+4]
    if strand=='-':
        ref7mer = str(Seq(ref7mer).reverse_complement())

    res = re.search(test_motif, ref7mer)
    if res is None:
        return None, None, None

    span = res.span()
    shift = span[0] - 1
    if strand=='+':
        chromStart += shift
        chromEnd += shift
    else:
        chromStart -= shift
        chromEnd -= shift

    ref5mer = this_ref[chrom][chromStart - 2:chromStart + 3]
    if strand=='-':
        ref5mer = str(Seq(ref5mer).reverse_complement())

    return chromStart, chromEnd, ref5mer


df_match = []
df_shifted = []
df_dumped = []
for _, row in df_in.iterrows():
    chrom = row['chr'].lstrip('chr')
    chromEnd = row['pos']
    chromStart = chromEnd - 1
    strand = row['strand']

    ref5mer = ref[chrom][chromStart-2:chromStart+3]
    if strand=='-':
        ref5mer = str(Seq(ref5mer).reverse_complement())

    if ref5mer==row['Motif_1']:
        df_match.append([chrom, chromStart, chromEnd, 'psi', round(row['Frac_Ave %'], 1), strand, ref5mer])
    else:
        chromStart, chromEnd, ref5mer = expanded_search(row, ref, row['Motif_1'])
        if (chromStart is None) or (ref5mer!=row['Motif_1']):
            df_dumped.append(row)
            continue
        else:
            df_shifted.append([chrom, chromStart, chromEnd, 'psi', round(row['Frac_Ave %'], 1), strand, ref5mer])

df_match = pd.DataFrame(df_match, columns=bed_fields)
df_shifted = pd.DataFrame(df_shifted, columns=bed_fields)
df_dumped = pd.DataFrame(df_dumped)

df_out = pd.concat((df_match, df_shifted))

df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False)