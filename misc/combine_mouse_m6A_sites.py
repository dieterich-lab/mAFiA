import os
import numpy as np
import pandas as pd
from Bio import SeqIO
col_names = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

ref_fasta = '/home/adrian/Data/genomes/GRCm38_102.fa'
ref = {}
for record in SeqIO.parse(ref_fasta, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = record.seq
    elif record.id == 'MT':
        ref['M'] = record.seq

bed_dir = '/home/adrian/Data/Dewenter_TAC_Backs_lab/achan/bed_files'

df_sham = pd.read_csv(
    os.path.join(bed_dir, 'nochr_predicted_m6A_sham.bed'),
    sep='\t',
    names=col_names,
    dtype={'chrom': str}
)
df_sham['name'] = 'm6A'
df_sham['score'] = np.round(df_sham['score'], 3)

df_tac = pd.read_csv(
    os.path.join(bed_dir, 'nochr_predicted_m6A_TAC.bed'),
    sep='\t',
    names=col_names,
    dtype={'chrom': str}
)
df_tac['name'] = 'm6A'
df_tac['score'] = np.round(df_tac['score'], 3)

df_merged = pd.merge(df_sham, df_tac, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand'], suffixes=['_sham', '_tac'])

ref5mer = []
for _, row in df_merged.iterrows():
    this_5mer = ref[row['chrom']][row['chromStart']-2:row['chromStart']+3]
    if row['strand'] == '-':
        this_5mer = this_5mer.reverse_complement()
    ref5mer.append(str(this_5mer))
df_merged['ref5mer'] = ref5mer

df_merged.to_csv(os.path.join(bed_dir, 'm6Anet_sites.bed'), index=False, sep='\t')