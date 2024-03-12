import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
from tqdm import tqdm
import pandas as pd

ref_file = '/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
# ref_file = '/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa'

out_dir = '/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102'
# out_dir = '/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102'

os.makedirs(out_dir, exist_ok=True)

ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    # if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
    if (record.id.isnumeric()) or (record.id in ['X']):
        ref[record.id] = str(record.seq)

all_chroms = [str(x) for x in ref.keys()]

# all_chroms = [str(xx) for xx in sorted([int(x) for x in ref.keys() if x.isnumeric()])]
# all_chroms.extend(['X', 'Y', 'MT'])

# DRACH_motifs = [
#     'AGACT',
#     'GAACT',
#     'GGACA',
#     'GGACC',
#     'GGACT',
#     'TGACT',
# ]

# mod_motifs = {
#     'm6A': [
#         'AAACA',
#         'AAACC',
#         'AAACT',
#         'AGACA',
#         'AGACC',
#         'AGACT',
#         'GAACA',
#         'GAACC',
#         'GAACT',
#         'GGACA',
#         'GGACC',
#         'GGACT',
#         'TAACA',
#         'TAACC',
#         'TAACT',
#         'TGACA',
#         'TGACC',
#         'TGACT'
#     ],
#     'psi': [
#         'AGTGG',
#         'GGTCC',
#         'GGTGG',
#         'GTTCA',
#         'GTTCC',
#         'GTTCG',
#         'TGTAG',
#         'TGTGG'
#     ]
# }

mod_motifs = {
    'Gm': [
        "ACGTC",
        "ATGAC",
        "ATGCT",
        "ATGTC",
        "ATGTG",
        "ATGTT",
        "CCGCC",
        "CTGCC",
        "CTGCG",
        "CTGCT",
        "CTGTA",
        "CTGTC",
        "CTGTG",
        "GCGCC",
        "GTGCA",
        "GTGCC",
        "GTGTC",
        "TCGCC"
    ],
}

# DRACH_motifs = [
#     'AAACA',
#     'AAACC',
#     'AAACT',
#     'AGACA',
#     'AGACC',
#     'AGACT',
#     'GAACA',
#     'GAACC',
#     'GAACT',
#     'GGACA',
#     'GGACC',
#     'GGACT',
#     'TAACA',
#     'TAACC',
#     'TAACT',
#     'TGACA',
#     'TGACC',
#     'TGACT',
# ]

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
    chromStart_name_strand_ref5mer = []
    for this_mod, this_mod_motifs in mod_motifs.items():
        for this_motif in this_mod_motifs:
            for this_match in re.finditer(this_motif, seq_forward):
                chromStart_name_strand_ref5mer.append((this_match.start()+2, this_mod, '+', this_motif))
            for this_match in re.finditer(this_motif, seq_revcomp):
                chromStart_name_strand_ref5mer.append((seq_len-this_match.start()-3, this_mod, '-', this_motif))
    chromStart_name_strand_ref5mer.sort()

    this_chrom_df = pd.DataFrame(chromStart_name_strand_ref5mer, columns=['chromStart', 'name', 'strand', 'ref5mer'])
    this_chrom_df['chrom'] = this_chrom
    this_chrom_df['chromEnd'] = this_chrom_df['chromStart'] + 1
    this_chrom_df['score'] = '.'
    this_chrom_df = this_chrom_df[fields]

    print(f'Writing out {len(this_chrom_df)} sites')
    this_chrom_df.to_csv(os.path.join(out_dir, f'Gm.GRCm38_102.chr{this_chrom}.bed'), sep='\t', index=False, header=True)