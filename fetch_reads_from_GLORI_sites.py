import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/GLORI'))
import argparse
from tqdm import tqdm
from glob import glob
import pandas as pd
import torch
from models import objectview
import pysam
from Bio import SeqIO
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import load_model, get_nucleotides_aligned_to_site
from feature_classifiers import get_outlier_ratio_from_features
from time import time

parser = argparse.ArgumentParser()
parser.add_argument('--perc_thresh', help='percolation threshold')
args = parser.parse_args()
PERC_THRESH = float(args.perc_thresh)
print('Clustering at threshold {:.2f}'.format(PERC_THRESH))

glori_file = os.path.join(HOME, 'Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
df_glori = pd.read_csv(glori_file)

ref_file = os.path.join(HOME, 'Data/genomes/GRCh38_96.fa')
ref = {}
print('Parsing genome...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

### WT ###
wt_bam_file = os.path.join(HOME, 'inference/HEK293_IVT_2_q50_10M/epoch29_HEK293A_WT_aligned_to_genome.bam.sorted')
wt_bam = pysam.AlignmentFile(wt_bam_file, 'rb')
wt_fast5_dir = '/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all'
wt_f5_paths = glob(os.path.join(wt_fast5_dir, '*.fast5'), recursive=True)
wt_index_read_ids = {}
print('Parsing WT fast5 files...', flush=True)
for f5_filepath in wt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        wt_index_read_ids[read_id] = f5_filepath
print('{} WT reads collected'.format(len(wt_index_read_ids)), flush=True)

### IVT ###
ivt_bam_file = os.path.join(HOME, 'inference/HEK293_IVT_2_q50_10M/epoch29_HEK293_IVT_2_aligned_to_genome.bam.sorted')
ivt_bam = pysam.AlignmentFile(ivt_bam_file, 'rb')
ivt_fast5_dir = '/home/achan/Data/Isabel_IVT_Nanopore/HEK293_IVT_2/fast5_pass'
ivt_f5_paths = glob(os.path.join(ivt_fast5_dir, '*.fast5'), recursive=True)
ivt_index_read_ids = {}
print('Parsing IVT fast5 files...', flush=True)
for f5_filepath in ivt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        ivt_index_read_ids[read_id] = f5_filepath
print('{} IVT reads collected'.format(len(ivt_index_read_ids)), flush=True)

### load model, device ###
model_path = os.path.join(HOME, 'pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch')
torchdict = torch.load(model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(model_path, fixed_config)

### search by GLORI sites ###
MIN_COVERAGE = 50
# PERC_THRESH = 0.7
outfile = os.path.join(HOME, 'Data/GLORI/df_outlier_ratios_thresh{:.2f}.tsv'.format(PERC_THRESH))
# print('Going through GLORI m6A sites...', flush=True)
df_outlier = pd.DataFrame()
counts = 0
for ind, row in df_glori.iterrows():
    ### debug ###
    # ind = 6047
    # row = df_glori.iloc[ind]
    # ind = 1466
    # row = df_glori.iloc[ind]
    # ind = 1715
    # row = df_glori.iloc[ind]
    # ind = 2299
    # row = df_glori.iloc[ind]
    # ind = 4112
    # row = df_glori.iloc[ind]
    # ind = 2298
    # row = df_glori.iloc[ind]

    print('\nSite {}'.format(ind), flush=True)
    chr = row['Chr'].lstrip('chr')
    strand = row['Strand']
    site = row['Sites'] - 1   # 0-based
    glori_ratio = row['Ratio']
    ref_motif = ref[chr][site-2:site+3]

    if not (((chr.isnumeric()) or (chr in ['X', 'Y'])) and (strand=='+')):
        continue

    # print('Collecting WT features...', flush=True)
    # tic = time()
    # v1_wt_site_motif_features = collect_features_from_aligned_site(fixed_model, fixed_device, fixed_config, wt_bam, wt_index_read_ids, chr, site, MIN_COVERAGE)
    # elapsed = time() - tic
    # print('v1 time elapsed: {:.1f} secs'.format(elapsed))

    # tic = time()
    wt_site_motif_features = get_nucleotides_aligned_to_site(fixed_model, fixed_device, fixed_config, wt_bam, wt_index_read_ids, chr, site, MIN_COVERAGE)
    # elapsed = time() - tic
    # print('v2 time elapsed: {:.1f} secs'.format(elapsed))

    # print('Collecting IVT features...', flush=True)
    ivt_site_motif_features = get_nucleotides_aligned_to_site(fixed_model, fixed_device, fixed_config, ivt_bam, ivt_index_read_ids, chr, site, MIN_COVERAGE, enforce_motif=ref_motif)

    if (len(wt_site_motif_features)>MIN_COVERAGE) and (len(ivt_site_motif_features)>MIN_COVERAGE):
        print('=========================================================', flush=True)
        print('chr{}, pos{}, strand{}'.format(chr, site, strand), flush=True)
        print('Reference motif {}'.format(ref_motif), flush=True)
        # print('Mod. ratio = {:.2f}'.format(glori_ratio), flush=True)
        print('{} feature vectors collected from WT'.format(len(wt_site_motif_features)), flush=True)
        print('{} feature vectors collected from IVT'.format(len(ivt_site_motif_features)), flush=True)
        print('Now clustering features...', flush=True)
        outlier_ratio = get_outlier_ratio_from_features(ivt_site_motif_features, wt_site_motif_features, ref_motif, PERC_THRESH)
        print('Calculated outlier ratio {:.2f} [GLORI {:.2f}]'.format(outlier_ratio, glori_ratio), flush=True)
        print('=========================================================', flush=True)
        if outlier_ratio!=-1:
            new_row = row.copy()
            new_row['motif'] = ref_motif
            new_row['ratio_outlier'] = outlier_ratio
            new_row['num_features_ivt'] = len(ivt_site_motif_features)
            new_row['num_features_wt'] = len(wt_site_motif_features)
            df_outlier = pd.concat([df_outlier, new_row.to_frame().T])
            counts += 1
            if counts%5==0:
                df_outlier.to_csv(outfile, sep='\t')
