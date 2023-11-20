from Bio import SeqIO
import re
import pandas as pd

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

all_chroms = [str(xx) for xx in sorted([int(x) for x in ref.keys() if x.isnumeric()])]
all_chroms.extend(['X', 'Y', 'MT'])

DRACH_motifs = [
    'AAACA',
    'AAACC',
    'AAACT',
    'AGACA',
    'AGACC',
    'AGACT',
    'GAACA',
    'GAACC',
    'GAACT',
    'GGACA',
    'GGACC',
    'GGACT',
    'TAACA',
    'TAACC',
    'TAACT',
    'TGACA',
    'TGACC',
    'TGACT',
]

for this_chrom in all_chroms:
