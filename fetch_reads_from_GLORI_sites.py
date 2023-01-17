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
    # print(record.id)
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

bam_file = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype/minimap/aligned_genome.bam.sorted')
bam = pysam.AlignmentFile(bam_file, 'rb')

# motifs = []
for _, row in df_glori.iterrows():
    chr = row['Chr'].lstrip('chr')

    if not ((chr.isnumeric()) or (chr in ['X', 'Y'])):
        continue

    strand = row['Strand']
    site = row['Sites'] - 1   # 0-based
    ratio = row['Ratio']

    if strand=='-':
        continue

    ref_motif = ref[chr][site-2:site+3]
    if strand=='-':
        ref_motif = str(Seq(ref_motif).reverse_complement())

    # overlap_reads = bam.fetch(chr, site, site+1)
    # for this_read in overlap_reads:
    #     print(str(this_read))

    for pileupcolumn in bam.pileup(chr, site, site+1, truncate=True, min_base_quality=0):
        if pileupcolumn.pos == site:
            # coverage = pileupcolumn.n
            coverage = pileupcolumn.get_num_aligned()
            if coverage>100:
                print('\nchr {}, pos {}, strand {}'.format(chr, pileupcolumn.pos, strand))
                print('Reference motif = {}'.format(ref_motif))
                print('Mod. ratio = {}'.format(ratio))
                print('coverage = {}'.format(coverage))
                valid_counts = 0
                for pileupread in pileupcolumn.pileups:
                    # query_position = pileupread.query_position_or_next
                    query_position = pileupread.query_position
                    if query_position and pileupread.alignment.query_sequence[query_position]=='A':
                        valid_counts += 1
                        print('\tmotif in read {} = {}, pos {}'.format(
                            pileupread.alignment.query_name,
                            pileupread.alignment.query_sequence[query_position-2:query_position+3],
                            query_position
                        ))
                print('Valid reads = {}'.format(valid_counts))

    # if strand=='+':
    #     # print(ref_seq)
    #     motifs.append(ref_seq)
    # else:
    #     # print(Seq(ref_seq).reverse_complement())
    #     motifs.append(str(Seq(ref_seq).reverse_complement()))