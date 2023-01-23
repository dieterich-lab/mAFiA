import os
HOME = os.path.expanduser('~')
from glob import glob
from utils import med_mad
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

bam_file = os.path.join(HOME, 'inference/HEK293_IVT_2_q50_10M/epoch29_HEK293A_WT_aligned_to_genome.bam.sorted')
bam = pysam.AlignmentFile(bam_file, 'rb')

# query_read_ids = []
# for read in bam.fetch():
#     if read.flag==0:
#         query_read_ids.append(read.query_name)

### index read ids ###
print('Indexing fast5 reads...')
fast5_dir = '/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all'
f5_paths = glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)

index_read_ids = {}
for f5_filepath in f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        index_read_ids[read_id] = f5_filepath

### parallel version ###
# from threading import Thread
# from multiprocessing import Process
#
# def worker(f5_filepath, wanted_read_ids, id_signal):
#     try:
#         print('Reading {}'.format(f5_filepath))
#         f5 = get_fast5_file(f5_filepath, mode="r")
#         for read_id in tqdm(f5.get_read_ids()):
#             # index_read_ids[read_id] = f5_filepath
#             if read_id in wanted_read_ids:
#                 read = f5.get_read(read_id)
#                 signal = read.get_raw_data(scale=True)
#                 med, mad = med_mad(signal)
#                 id_signal[read_id] = (signal - med) / mad
#     except:
#         print('Error with {}'.format(f5_filepath))
#     return True
#
# id_signal = [{} for x in f5_paths]
# processes = []
# for ii, f5_filepath in enumerate(f5_paths):
#     # process = Thread(target=worker, args=[f5_filepath, query_read_ids, id_signal[ii]])
#     p = Process(target=worker, args=[f5_filepath, query_read_ids, id_signal[ii]])
#     p.start()
#     processes.append(p)
#
# for p in processes:
#     p.join()
#
# id_signal = {k: v for d in id_signal for (k, v) in d.items()}

### serial version ###
# for f5_filepath in tqdm(glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)):
#     f5 = get_fast5_file(f5_filepath, mode="r")
#     for read_id in f5.get_read_ids():
#         # index_read_ids[read_id] = f5_filepath
#         if read_id in query_read_ids:
#             read = f5.get_read(read_id)
#             signal = read.get_raw_data(scale=True)
#             med, mad = med_mad(signal)
#             id_signal[read_id] = (signal - med) / mad

### search by GLORI sites ###
print('Going through GLORI m6A sites...')
for ind, row in df_glori.iterrows():
    ### debug ###
    ind = 6047
    row = df_glori.iloc[ind]

    chr = row['Chr'].lstrip('chr')
    strand = row['Strand']
    site = row['Sites'] - 1   # 0-based
    glori_ratio = row['Ratio']

    if not ((chr.isnumeric()) or (chr in ['X', 'Y']) or (strand=='+')):
        continue

    ref_motif = ref[chr][site-2:site+3]

    site_features = {}
    for pileupcolumn in bam.pileup(chr, site, site+1, truncate=True, min_base_quality=0):
        if pileupcolumn.pos == site:
            coverage = pileupcolumn.get_num_aligned()
            if coverage>100:
                valid_counts = 0
                print('\nCollecting features for site {}, chr{}, pos{}, strand{}'.format(ind, chr, pileupcolumn.pos, strand))
                for pileupread in tqdm(pileupcolumn.pileups):
                    query_name = pileupread.alignment.query_name
                    query_position = pileupread.query_position_or_next
                    # query_position = pileupread.query_position
                    flag = pileupread.alignment.flag

                    # if query_position and (flag==0) and (query_name in index_read_ids.keys()) and (pileupread.alignment.query_sequence[query_position] == 'A'):
                    if query_position and (flag==0) and (query_name in index_read_ids.keys()) and (glori_ratio>0.9):
                    # if query_position and (flag==0) and (query_name in id_signal.keys()):
                        valid_counts += 1
                        query_motif = pileupread.alignment.query_sequence[query_position-2:query_position+3]
                        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
                        # this_read_signal = id_signal[query_name]
                        this_read_features, the_read_motif = extract_features_from_signal(this_read_signal, query_position, query_motif)
                        if this_read_features is not None:
                            site_features[query_name] = this_read_features

                if valid_counts>0:
                    print('Reference motif = {}'.format(ref_motif))
                    print('Mod. ratio = {}'.format(glori_ratio))
                    print('coverage = {}'.format(coverage))
                    # print('Valid reads = {}'.format(valid_counts))
                    print('{} feature vectors collected'.format(len(site_features)))
                else:
                    print('No valid reads!')

    if len(site_features)>0:
        print('Now clustering features...')
        outlier_ratio = cluster_features(site_features)
        print('Calculated outlier ratio {:.3f} [GLORI {:.3f}]'.format(outlier_ratio, glori_ratio))