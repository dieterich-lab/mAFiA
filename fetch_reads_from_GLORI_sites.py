import os
HOME = os.path.expanduser('~')
from glob import glob
from tqdm import tqdm
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from ont_fast5_api.fast5_interface import get_fast5_file
from utils import get_norm_signal_from_read_id
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

# bam_file = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype/minimap/aligned_genome.bam.sorted')
bam_file = os.path.join(HOME, 'inference/HEK293_IVT_2_q50_10M/epoch29_HEK293A_WT_aligned_to_genome.bam.sorted')
bam = pysam.AlignmentFile(bam_file, 'rb')

### index read ids ###
# fast5_dir = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype/fast5_filtered')
fast5_dir = '/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all'
# id_signal = {}
index_read_ids = {}
print('Indexing read ids...')
for f5_filepath in tqdm(glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)):
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        index_read_ids[read_id] = f5_filepath

### search by GLORI sites ###
for ind, row in tqdm(df_glori.iterrows()):
    # ### debug ###
    # ind = 2746
    # row = df_glori.iloc[ind]

    chr = row['Chr'].lstrip('chr')
    strand = row['Strand']
    site = row['Sites'] - 1   # 0-based
    ratio = row['Ratio']

    # if not ((chr.isnumeric()) or (chr in ['X', 'Y']) or (strand=='+')):
    #     continue

    ref_motif = ref[chr][site-2:site+3]
    # if strand=='-':
    #     ref_motif = str(Seq(ref_motif).reverse_complement())

    # overlap_reads = bam.fetch(chr, site, site+1)
    # for this_read in overlap_reads:
    #     print(str(this_read))

    site_features = {}
    for pileupcolumn in bam.pileup(chr, site, site+1, truncate=True, min_base_quality=0):
        if pileupcolumn.pos == site:
            # coverage = pileupcolumn.n
            coverage = pileupcolumn.get_num_aligned()
            if coverage>100:
                valid_counts = 0
                for pileupread in pileupcolumn.pileups:
                    query_name = pileupread.alignment.query_name
                    # query_position = pileupread.query_position_or_next
                    query_position = pileupread.query_position
                    flag = pileupread.alignment.flag
                    # print(query_name, query_position, flag, pileupread.alignment.query_sequence[query_position-2:query_position+3])
                    # if (flag!=1):
                    #     continue
                    # if flag==16:
                    #     reverse_complement = True
                    # else:
                    #     reverse_complement = False
                    if query_position and (flag==0) and (query_name in index_read_ids.keys()) and (pileupread.alignment.query_sequence[query_position] == 'A'):
                    # if query_position and (query_name in index_read_ids.keys()) and (pileupread.alignment.query_sequence[query_position] == 'A'):
                        valid_counts += 1
                        # print('\tmotif in read {} = {}, pos {}'.format(
                        #     query_name,
                        #     pileupread.alignment.query_sequence[query_position-2:query_position+3],
                        #     query_position
                        # ))
                        query_motif = pileupread.alignment.query_sequence[query_position-2:query_position+3]
                        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
                        this_read_features = extract_features_from_signal(this_read_signal, query_position, query_motif)
                        # print('{}, pos {}'.format(query_name, query_position))
                        # print(this_read_features)
                        if this_read_features is not None:
                            site_features[query_name] = this_read_features
                        # else:
                        #     print('Query read not in fast5 directory!')

                if valid_counts>0:
                    print('\nrow {}, chr {}, pos {}, strand {}'.format(ind, chr, pileupcolumn.pos, strand))
                    print('Reference motif = {}'.format(ref_motif))
                    print('Mod. ratio = {}'.format(ratio))
                    print('coverage = {}'.format(coverage))
                    print('Valid reads = {}'.format(valid_counts))
                    print('{} features collected'.format(len(site_features)))

    # outlier_ratio = cluster_features(collected_features)

    # if strand=='+':
    #     # print(ref_seq)
    #     motifs.append(ref_seq)
    # else:
    #     # print(Seq(ref_seq).reverse_complement())
    #     motifs.append(str(Seq(ref_seq).reverse_complement()))