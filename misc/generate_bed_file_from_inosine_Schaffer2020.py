import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
import numpy as np


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

in_file = '/home/adrian/Data/inosine/Schaffer2020_Supplementary_Tables_S2_HEK293_GRCh38.tsv'
df_in = pd.read_csv(in_file, sep='\t', dtype={'chrom': str})
out_file = '/home/adrian/Data/inosine/HEK293_Schaffer2020.bed'

df_in['chrom'] = [this_chrom.lstrip('chr') for this_chrom in df_in['chrom']]
df_in = df_in[
    (df_in['chrom'].isin(ref.keys()))
    * ~(df_in['HEK293_a-%editing'].isna() & df_in['HEK293_a-%editing'].isna())
]

df_out = []
for _, row in df_in.iterrows():
    chrom = row['chrom']
    chromStart = row['chromStart']
    chromEnd = row['chromEnd']

    editing_ratios = row[['HEK293_a-%editing', 'HEK293_b-%editing']].values
    modRatio = np.mean(editing_ratios[~np.isnan(np.float64(editing_ratios))]) * 100.0
    # modRatio = (row['HEK293_a-%editing'] + row['HEK293_b-%editing']) / 2 * 100.0

    ref5mer = ref[chrom][chromStart - 2:chromStart + 3]

    if ref5mer[2] == 'A':
        strand = '+'
    elif ref5mer[2] == 'T':
        strand = '-'
        ref5mer = str(Seq(ref5mer).reverse_complement())
    else:
        continue

    print(strand, ref5mer)
    df_out.append([chrom, chromStart, chromEnd, 'ino', modRatio, strand, ref5mer])
df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False, float_format='%.1f')