import os
HOME = os.path.expanduser('~')
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

glori_file = os.path.join(HOME, 'Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
df_glori = pd.read_csv(glori_file)

ref_file = os.path.join(HOME, 'Data/genomes/GRCh38_96.fa')
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    print(record.id)
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

bam_file = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype/minimap/aligned_genome.bam.sorted')
bam = pysam.AlignmentFile(bam_file, 'rb')

motifs = []
for _, row in df_glori.iterrows():
    chr = row['Chr'].lstrip('chr')

    if not ((chr.isnumeric()) or (chr in ['X', 'Y'])):
        continue

    strand = row['Strand']
    site = row['Sites'] - 1   # 0-based
    ratio = row['Ratio']

    ref_seq = ref[chr][site-2:site+3]
    if strand=='+':
        # print(ref_seq)
        motifs.append(ref_seq)
    else:
        # print(Seq(ref_seq).reverse_complement())
        motifs.append(str(Seq(ref_seq).reverse_complement()))

    # overlap_reads = bam.fetch(chr, site, site+1)
    # for this_read in overlap_reads:
    #     print(str(this_read))
