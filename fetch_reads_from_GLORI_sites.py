import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/GLORI'))
from tqdm import tqdm
from glob import glob
import pandas as pd
import pysam
from Bio import SeqIO
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import collect_features_from_aligned_site
from cluster_features import get_outlier_ratio_from_features

glori_file = os.path.join(HOME, 'Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
df_glori = pd.read_csv(glori_file)

ref_file = os.path.join(HOME, 'Data/genomes/GRCh38_96.fa')
ref = {}
print('Parsing genome...')
for record in tqdm(SeqIO.parse(ref_file, 'fasta')):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

### WT ###
wt_bam_file = os.path.join(HOME, 'inference/HEK293_IVT_2_q50_10M/epoch29_HEK293A_WT_aligned_to_genome.bam.sorted')
wt_bam = pysam.AlignmentFile(wt_bam_file, 'rb')
wt_fast5_dir = '/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all'
wt_f5_paths = glob(os.path.join(wt_fast5_dir, '*.fast5'), recursive=True)
wt_index_read_ids = {}
print('Parsing WT fast5 files...')
for f5_filepath in tqdm(wt_f5_paths):
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        wt_index_read_ids[read_id] = f5_filepath
print('{} WT reads collected'.format(len(wt_index_read_ids)))

### IVT ###
ivt_bam_file = os.path.join(HOME, 'inference/HEK293_IVT_2_q50_10M/epoch29_HEK293_IVT_2_aligned_to_genome.bam.sorted')
ivt_bam = pysam.AlignmentFile(ivt_bam_file, 'rb')
ivt_fast5_dir = '/home/achan/Data/Isabel_IVT_Nanopore/HEK293_IVT_2/fast5_pass'
ivt_f5_paths = glob(os.path.join(ivt_fast5_dir, '*.fast5'), recursive=True)
ivt_index_read_ids = {}
print('Parsing IVT fast5 files...')
for f5_filepath in tqdm(ivt_f5_paths):
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        ivt_index_read_ids[read_id] = f5_filepath
print('{} IVT reads collected'.format(len(ivt_index_read_ids)))

### search by GLORI sites ###
MIN_COVERAGE = 50
PERC_THRESH = 0.85
# print('Going through GLORI m6A sites...')
for ind, row in df_glori.iterrows():
    ### debug ###
    # ind = 6047
    # row = df_glori.iloc[ind]
    ind = 1466
    row = df_glori.iloc[ind]

    chr = row['Chr'].lstrip('chr')
    strand = row['Strand']
    site = row['Sites'] - 1   # 0-based
    glori_ratio = row['Ratio']
    ref_motif = ref[chr][site-2:site+3]

    if not (((chr.isnumeric()) or (chr in ['X', 'Y'])) and (strand=='+')):
        continue

    # print('Collecting WT features...')
    wt_site_motif_features = collect_features_from_aligned_site(wt_bam, wt_index_read_ids, chr, site, MIN_COVERAGE)
    # print('Collecting IVT features...')
    ivt_site_motif_features = collect_features_from_aligned_site(ivt_bam, ivt_index_read_ids, chr, site, MIN_COVERAGE)

    if (len(wt_site_motif_features)>MIN_COVERAGE) and (len(ivt_site_motif_features)>MIN_COVERAGE):
        print('\nSite {}, chr{}, pos{}, strand{}'.format(ind, chr, site, strand))
        print('Reference motif {}'.format(ref_motif))
        # print('Mod. ratio = {:.2f}'.format(glori_ratio))
        print('{} feature vectors collected from WT'.format(len(wt_site_motif_features)))
        print('{} feature vectors collected from IVT'.format(len(ivt_site_motif_features)))
        print('Now clustering features...')
        outlier_ratio = get_outlier_ratio_from_features(ivt_site_motif_features, wt_site_motif_features, ref_motif, PERC_THRESH)
        print('Calculated outlier ratio {:.2f} [GLORI {:.2f}]'.format(outlier_ratio, glori_ratio))