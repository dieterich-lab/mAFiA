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

in_file = '/home/adrian/Data/m5C-TAC-Seq/mmc5.xlsx'
df_in = pd.read_excel(in_file,
                      sheet_name='HEK293T_cutoff-1',
                      usecols=[
                          'Chromosome',
                          'Position',
                          'Strand',
                          'C-to-T_ratio_m5C-TAC_rep1 (%)',
                          'C-to-T_ratio_m5C-TAC_rep2 (%)'
                      ])
out_file = '/home/adrian/Data/m5C-TAC-Seq/HEK293T_cutoff-1.bed'

df_in['chrom'] = ['MT' if this_chrom == 'chrM' else this_chrom.lstrip('chr') for this_chrom in df_in['Chromosome']]

df_match = []
for _, row in df_in.iterrows():
    chrom = row['chrom']
    chromEnd = row['Position']
    chromStart = chromEnd - 1
    strand = row['Strand']

    ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    if strand == '-':
        ref5mer = str(Seq(ref5mer).reverse_complement())
    print(strand, ref5mer)

    test_chars = ref5mer[2], ref5mer[1], ref5mer[3]
    if [test_chars[i]=='T' for i in range(span)]:
        chromEnd = chromStart + 1
        avg_deletion_ration = round((row['rep1_deletion_ratio'] + row['rep2_deletion_ratio']) * 0.5 * 100.0, 1)
        df_match.append([chrom, chromStart, chromEnd, 'psi', avg_deletion_ration, strand, ref5mer])


df_match = pd.DataFrame(df_match, columns=bed_fields)
# df_consensus = pd.DataFrame(df_consensus, columns=bed_fields)
# df_out = pd.concat([df_match, df_consensus])
df_out = df_match
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False, float_format='%.1f')