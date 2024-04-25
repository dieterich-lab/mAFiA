import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
from liftover import get_lifter

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
]

in_file = '/home/adrian/Data/Psi-Seq/mmc3.xlsx'
ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
sheet_name = 'AllHumanSites'
out_file = '/home/adrian/Data/Psi-Seq/Psi-Seq_HEK293.bed'

ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)


df_in = pd.read_excel(in_file, sheet_name=sheet_name,
                      usecols=['gcoords', 'strand', 'SequenceSurroundingSite', 'siControlrep1_PSIratio', 'siControlrep2_PSIratio'])


df_out = []
converter = get_lifter('hg19', 'hg38', one_based=False)
for _, row in df_in.iterrows():
    gcoords_old = row['gcoords'].split(':')
    gcoords = converter[gcoords_old[0]][int(gcoords_old[1])][0][:2]
    chrom = gcoords[0].lstrip('chr')
    strand = row['strand']

    if strand=='+':
        chromEnd = int(gcoords[1]) - 1
        chromStart = chromEnd - 1
        refseq = ref[chrom][chromStart - 9:chromStart + 12]
    else:
        chromStart = int(gcoords[1])
        chromEnd = chromStart + 1
        refseq = ref[chrom][chromStart - 11:chromStart + 10]

    ref5mer = ref[chrom][chromStart-2:chromStart+3]
    if strand=='-':
        refseq = str(Seq(refseq).reverse_complement())
        ref5mer = str(Seq(ref5mer).reverse_complement())
    if refseq!=row['SequenceSurroundingSite']:
        print(strand, refseq, row['SequenceSurroundingSite'])
        continue

    diff_rep = row['siControlrep1_PSIratio'] - row['siControlrep2_PSIratio']
    # if diff_rep>=0.1:
    #     continue
    avg_rep = round(0.5 * (row['siControlrep1_PSIratio'] + row['siControlrep2_PSIratio']) * 100.0, 1)

    df_out.append([chrom, chromStart, chromEnd, 'psi', avg_rep, strand, ref5mer])
df_out = pd.DataFrame(df_out, columns=bed_fields)

df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False)