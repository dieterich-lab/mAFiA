import pandas as pd
pd.set_option('display.max_columns', 500)
from Bio import SeqIO
from Bio.Seq import Seq
import re
import numpy as np

input_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Nm-Mut-seq_Supp_Tables.xlsx'
# sheet_name = 'S1_HeLa_known sites_rRNA_WT'
# sheet_name = 'S2_HeLa_full set rRNA_WT'
sheet_name = 'S6_HepG2_mRNA_WT'
df_in = pd.read_excel(input_file, sheet_name=sheet_name, skiprows=[0, 1])

# ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/rRNA_18S_28S.fasta'
# ref = {}
# with open(ref_file, 'r') as h_ref:
#     for record in SeqIO.parse(h_ref, 'fasta'):
#         ref[record.id] = str(record.seq)

motifs = [
    'AAGAA',
    'AAGAG',
    'AAGAT',
    'AAGCA',
    'AGGAA',
    'AGGAG',
    'AGGCC',
    'ATGGA',
    'CAGAA',
    'CAGAG',
    'CAGCA',
    'CAGCC',
    'CAGCT',
    'CAGGA',
    'CAGGC',
    'CAGGT',
    'CTGAA',
    'CTGAG',
    'CTGCA',
    'CTGCC',
    'CTGCT',
    'CTGGA',
    'CTGTG',
    'GAGAA',
    'GAGAG',
    'GAGCA',
    'GAGCC',
    'GAGCT',
    'GAGGA',
    'GAGGC',
    'GTGGA',
    'GTGGC',
    'TAGAG',
    'TAGCG',
    'TGGAA',
    'TGGAG',
    'TGGCA',
    'TGGCC',
    'TGGCT',
    'TGGTG',
    'TTGAA',
    'TTGGA'
]

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref = {}
with open(ref_file, 'r') as h_ref:
    for record in SeqIO.parse(h_ref, format='fasta'):
        if record.id.isnumeric() or record.id=='X':
            ref[record.id] = str(record.seq)

out_bed_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/HepG2_mRNA_WT_Gm_sites.bed'

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
    # 'origPosition'
]

# df_gm = df_in[df_in['Mod_Base']=='Gm']
df_gm = df_in[df_in['class']=='Gm']

df_out = []
# counts = 0
for _, row in df_gm.iterrows():
    chrom = row['chr'].lstrip('chr')
    # if chrom=='5.8S':
    #     continue

    chromStart = row['position'] - 1
    # name = row['Mod_Base']
    name = row['class']
    strand = row['strand']
    motif = row['motif']

    ref5mer = ref[chrom][(chromStart-2):(chromStart+3)]
    if strand=='-':
        ref5mer = str(Seq(ref5mer).reverse_complement())

    # if chrom=='28S':
    #     all_pos = []
    #     for this_find in re.finditer(ref5mer, ref[chrom]):
    #         all_pos.append(this_find.span()[0]+2)
    #     all_pos = np.array(all_pos)
    #     diff_vec = all_pos - chromStart
    #     this_diff = np.min(diff_vec[diff_vec>=0])
    #     if this_diff in [13, 21, 30]:
    #         chromStart = chromStart+this_diff
    # print(chrom, chromStart, ref5mer, ref[chrom][(chromStart-2):(chromStart+3)])

    chromEnd = chromStart + 1

    # print(motif, ref5mer)

    if ref5mer==motif:
        score = round(row['Frac_Ave %'], 1)
        df_out.append([chrom, chromStart, chromEnd, name, score, strand, ref5mer])

df_out = pd.DataFrame(df_out, columns=bed_fields)
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_bed_file, sep='\t', index=False)