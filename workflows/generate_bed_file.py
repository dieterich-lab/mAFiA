import pandas as pd
pd.set_option('display.max_columns', 500)
from Bio import SeqIO
import re
import numpy as np

input_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Nm-Mut-seq_Supp_Tables.xlsx'
df_in = pd.read_excel(input_file, sheet_name='S1_HeLa_known sites_rRNA_WT', skiprows=[0, 1])

ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/rRNA_18S_28S.fasta'
ref = {}
with open(ref_file, 'r') as h_ref:
    for record in SeqIO.parse(h_ref, 'fasta'):
        ref[record.id] = str(record.seq)

out_bed_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/HeLa_rRNA_Gm_sites.bed'

# ref['28S'] = ref['28S']

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    # 'coverage',
    # 'modRatio',
    'ref5mer',
    'origPosition'
]

df_gm = df_in[df_in['Mod_Base']=='Gm']

df_out = []
for _, row in df_gm.iterrows():
    chrom = row['chr']
    if chrom=='5.8S':
        continue

    chromStart = row['position'] - 1
    name = row['Mod_Base']
    ref5mer = row['motif']

    if chrom=='28S':
        all_pos = []
        for this_find in re.finditer(ref5mer, ref[chrom]):
            all_pos.append(this_find.span()[0]+2)
        all_pos = np.array(all_pos)
        diff_vec = all_pos - chromStart
        this_diff = np.min(diff_vec[diff_vec>=0])
        if this_diff in [13, 21, 30]:
            chromStart = chromStart+this_diff
    # print(chrom, chromStart, ref5mer, ref[chrom][(chromStart-2):(chromStart+3)])

    chromEnd = chromStart + 1

    if ref5mer==ref[chrom][(chromStart-2):(chromStart+3)]:
        df_out.append([chrom, chromStart, chromEnd, name, '.', '+', ref5mer, row['position']])

df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out.to_csv(out_bed_file, sep='\t', index=False)