import os
import pandas as pd
import numpy as np

glori_file = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv'
out_dir = '/home/adrian/Data/GLORI/bed_files'
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

dict_chr = {
    'M' : 'MT'
}

df = pd.read_csv(glori_file)
unique_chr = df['Chr'].unique()
unique_chr.sort()

for this_chr in unique_chr:
    df_chr = df[df['Chr']==this_chr]
    num_sites = len(df_chr)
    df_bed = pd.DataFrame()
    chr_num =  this_chr.lstrip('chr')
    if chr_num=='M': chr_num = 'MT'
    df_bed[bed_fields[0]] = [chr_num] * num_sites
    ends = df_chr['Sites'].values
    df_bed[bed_fields[1]] = ends - 1
    df_bed[bed_fields[2]] = ends
    df_bed[bed_fields[3]] = 'm6A'
    df_bed[bed_fields[4]] = np.round(df_chr['NormeRatio'].values * 100.0, 1)
    df_bed[bed_fields[5]] = df_chr['Strand'].values

    df_bed.to_csv(os.path.join(out_dir, f'GLORI.{this_chr}.tsv'), sep='\t', index=False)