import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
from tqdm import tqdm
import pandas as pd

# ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref_file = '/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

all_chroms = [str(xx) for xx in sorted([int(x) for x in ref.keys() if x.isnumeric()])]
all_chroms.extend(['X', 'Y', 'MT'])

out_dir = '/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/site_annotations'
os.makedirs(out_dir, exist_ok=True)

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

fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

for this_chrom in all_chroms:
    print(f'Parsing chr{this_chrom}')
    seq_forward = ref[this_chrom]
    seq_revcomp = str(Seq(ref[this_chrom]).reverse_complement())
    seq_len = len(seq_forward)
    chromStart_strand_ref5mer = []
    for this_motif in tqdm(DRACH_motifs):
        for this_match in re.finditer(this_motif, seq_forward):
            chromStart_strand_ref5mer.append((this_match.start()+2, '+', this_motif))
        for this_match in re.finditer(this_motif, seq_revcomp):
            chromStart_strand_ref5mer.append((seq_len-this_match.start()-3, '-', this_motif))
    chromStart_strand_ref5mer.sort()

    this_chrom_df = pd.DataFrame(chromStart_strand_ref5mer, columns=['chromStart', 'strand', 'ref5mer'])
    this_chrom_df['chrom'] = this_chrom
    this_chrom_df['chromEnd'] = this_chrom_df['chromStart'] + 1
    this_chrom_df['name'] = 'DRACH'
    this_chrom_df['score'] = '.'
    this_chrom_df = this_chrom_df[fields]

    print(f'Writing out {len(this_chrom_df)} sites')
    this_chrom_df.to_csv(os.path.join(out_dir, f'DRACH.GRCh38_102.chr{this_chrom}.bed'), sep='\t', index=False, header=True)