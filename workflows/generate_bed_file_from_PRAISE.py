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

in_file = '/home/adrian/Data/BID_seq/41589_2023_1304_MOESM3_ESM.xlsx'
df_in = pd.read_excel(in_file, skiprows=range(2),
                      sheet_name='Supplementary Dataset 2',
                      usecols=[
                          'rep1_deletion_ratio',
                          'rep2_deletion_ratio',
                          'chr_site'])
out_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/psU/site_annotations/PRAISE.bed'

df_out = []
for _, row in df_in.iterrows():
    if (row['chr_site']=='NONE') or (len(row['chr_site'].split('_'))>2):
        continue
    chr = row['chr_site'].split('_')[0]
    site = row['chr_site'].split('_')[-1]
    chrom = chr.lstrip('chr')
    if '-' in site:
        continue
    chromEnd = int(site)
    chromStart = chromEnd - 1
    # strand = row['strand']

    # ref5mer = ref[chrom][chromStart-2:chromStart+3]
    # if strand=='-':
    #     ref5mer = str(Seq(ref5mer).reverse_complement())

    avg_deletion_ration = round((row['rep1_deletion_ratio'] + row['rep2_deletion_ratio']) * 0.5 * 100.0, 1)
    ref5mer = ref[chrom][chromStart-2:chromStart+3]

    if ref5mer[2]=='T':
        strand = '+'
    elif ref5mer[2]=='A':
        ref5mer = str(Seq(ref5mer).reverse_complement())
        strand = '-'
    else:
        continue

    # if strand=='-':
    #     ref9mer = str(Seq(ref9mer).reverse_complement())
    #
    # span = re.search(row['Motif_1'], ref9mer).span()
    # if span[0]!=2:
    #     continue
    #     # shift = span[0] - 2
    #     # if strand=='+':
    #     #     chromStart += shift
    #     #     chromEnd += shift
    #     # else:
    #     #     chromStart -= shift
    #     #     chromEnd -= shift
    # ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    # if strand == '-':
    #     ref5mer = str(Seq(ref5mer).reverse_complement())

    # print(ref5mer==row['Motif_1'])

    df_out.append([chrom, chromStart, chromEnd, 'psU', avg_deletion_ration, strand, ref5mer])
df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False)