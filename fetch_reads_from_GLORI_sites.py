import os
HOME = os.path.expanduser('~')
from glob import glob
from tqdm import tqdm
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from utils import med_mad
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import extract_features_from_signal
from cluster_features import cluster_features

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

fast5_dir = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype/fast5_filtered')
id_signal = {}
for f5_file in tqdm(glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)):
    f5 = get_fast5_file(f5_file, mode="r")
    for read in f5.get_reads():
        signal = read.get_raw_data(scale=True)
        signal_start = 0
        signal_end = len(signal)
        med, mad = med_mad(signal[signal_start:signal_end])
        signal = (signal[signal_start:signal_end] - med) / mad
        id_signal[read.get_read_id()] = signal

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

    collected_features = {}
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
                    query_name = pileupread.alignment.query_name
                    # query_position = pileupread.query_position_or_next
                    query_position = pileupread.query_position
                    if query_position and pileupread.alignment.query_sequence[query_position]=='A':
                        valid_counts += 1
                        # print('\tmotif in read {} = {}, pos {}'.format(
                        #     query_name,
                        #     pileupread.alignment.query_sequence[query_position-2:query_position+3],
                        #     query_position
                        # ))
                        if query_name in id_signal.keys():
                            query_motif = pileupread.alignment.query_sequence[query_position-2:query_position+3]
                            this_read_features = extract_features_from_signal(id_signal[query_name], query_position, query_motif)
                            print('{}, pos {}:'.format(query_name, query_position))
                            print(this_read_features)
                            collected_features[query_name] = this_read_features
                        # else:
                        #     print('Query read not in fast5 directory!')
                print('Valid reads = {}'.format(valid_counts))

    # outlier_ratio = cluster_features(collected_features)

    # if strand=='+':
    #     # print(ref_seq)
    #     motifs.append(ref_seq)
    # else:
    #     # print(Seq(ref_seq).reverse_complement())
    #     motifs.append(str(Seq(ref_seq).reverse_complement()))