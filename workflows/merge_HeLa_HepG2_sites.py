import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
from Bio import SeqIO
from Bio.Seq import Seq

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref = {}
with open(ref_file, 'r') as h_ref:
    for record in SeqIO.parse(h_ref, format='fasta'):
        if record.id.isnumeric() or record.id=='X':
            ref[record.id] = str(record.seq)

input_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Nm-Mut-seq_Supp_Tables.xlsx'
output_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/Gm_sites_mRNA_HeLa_HepG2_merged.bed'
img_out = '/home/adrian/img_out/Gmorah'

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
    'modRatio'
]

def import_data_sheet(sheet_name):
    df =  pd.read_excel(input_file, sheet_name=sheet_name, skiprows=[0, 1],
                        usecols=['class', 'chr', 'position', 'motif', 'strand', 'Frac_Ave %'])
    df_gm = df[df['class']=='Gm']
    df_gm = df_gm.rename(columns={'class': 'name', 'chr': 'chrom', 'position': 'chromEnd', 'motif': 'ref5mer', 'Frac_Ave %': 'modRatio'})
    df_gm['chrom'] = [this_chr.lstrip('chr') for this_chr in df_gm['chrom']]
    df_gm['chromStart'] = df_gm['chromEnd'] - 1
    df_gm['score'] = '.'
    return df_gm[bed_fields]

df_hela = import_data_sheet(sheet_name='S5_HeLa_mRNA_WT')
df_hepg2 = import_data_sheet(sheet_name='S6_HepG2_mRNA_WT')
df_merge = pd.merge(df_hela, df_hepg2, on=bed_fields[:-1], suffixes=['_HeLa', '_HepG2'])
df_merge.sort_values(by=['chrom', 'chromStart'], inplace=True)
df_merge.to_csv(output_file, sep='\t', index=False)

plt.figure(figsize=(5, 5))
plt.plot(df_merge['modRatio_HeLa'], df_merge['modRatio_HepG2'], 'o')
plt.xlim([-5, 105])
plt.ylim([-5, 105])
plt.xlabel('HeLa', fontsize=12)
# plt.xlabel('NmMutSeq')
plt.ylabel('HepG2', fontsize=12)
plt.title(f'{len(df_merge)} Gm sites on mRNA', fontsize=15)
plt.savefig(os.path.join(img_out, 'scatter_HeLa_vs_HepG2_mRNA.png'), bbox_inches='tight')

for _, row in df_merge.iterrows():
    chrom = row['chrom']
    chromStart = row['chromStart']
    strand = row['strand']
    ref5mer = row['ref5mer']

    check5mer = ref[chrom][chromStart-2:chromStart+3]
    if strand=='-':
        check5mer = str(Seq(check5mer).reverse_complement())
    print(ref5mer, check5mer, ref5mer==check5mer)

