import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
# pd.set_option('display.max_columns', 500)
from tqdm import tqdm

# in_bed = '/home/adrian/Data/TRR319_RMaP_BaseCalling/RNA004/dorado/dummy.bed'
# out_bed = '/home/adrian/Data/TRR319_RMaP_BaseCalling/RNA004/dorado/dummy.mAFiA.bed'
# reference = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_108/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

in_bed = '/prj/TRR319_RMaP_BaseCalling/RNA004/dorado/RNA004_HEK293_WT_RTA.bed'
out_bed = '/prj/TRR319_RMaP_BaseCalling/RNA004/dorado/RNA004_HEK293_WT_RTA.mAFiA.bed'
reference = '/biodb/genomes/homo_sapiens/GRCh38_108/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

thresh_coverage = 20

ref = {}
with gzip.open(reference, 'rt') as h_ref:
    for record in SeqIO.parse(h_ref, 'fasta'):
        # if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        if (record.id.isnumeric()) or (record.id in ['X']):
            ref[record.id] = str(record.seq)

cols_bedMethyl = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'thickStart',
    'thickEnd',
    'color',
    'coverage',
    'modRatio'
]
sel_ind = list(range(4)) + [5, 9, 10]
sel_cols = [cols_bedMethyl[i] for i in sel_ind]

cols_mAFiA = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
    'modRatio',
    'coverage',
]

dict_mods = {
    'a': 'm6A',
    '17802': 'psi'
}

df_in = pd.read_csv(in_bed, usecols=sel_ind, names=sel_cols, sep='\t', dtype={'chrom': str})
df_in_thresh = df_in[df_in['coverage']>=thresh_coverage]
for k, v in dict_mods.items():
    df_in_thresh.loc[df_in_thresh['name']==k, 'name'] = v

df_out = []
for _, row in tqdm(df_in_thresh.iterrows()):
    if not row['chrom'] in ref.keys():
        continue
    ref5mer = ref[row['chrom']][row['chromStart']-2:row['chromStart']+3]
    if row['strand']=='-':
        ref5mer = str(Seq(ref5mer).reverse_complement())
    row['ref5mer'] = ref5mer
    df_out.append(row)
df_out = pd.DataFrame(df_out)
df_out['score'] = '.'
df_out = df_out[cols_mAFiA]

df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_bed, sep='\t', header=True, index=False, float_format='%.1f')
