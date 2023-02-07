import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import argparse
from tqdm import tqdm
from glob import glob
import pandas as pd
import torch
from models import objectview
import pysam
from Bio import SeqIO
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import load_model, collect_features_from_aligned_site, collect_features_from_aligned_site_v2
from cluster_features import get_outlier_ratio_from_features
from time import time

# parser = argparse.ArgumentParser()
# parser.add_argument('--perc_thresh', help='percolation threshold')
# args = parser.parse_args()
# PERC_THRESH = float(args.perc_thresh)
# print('Clustering at threshold {:.2f}'.format(PERC_THRESH))
PERC_THRESH = 0.8

mod_file = os.path.join(HOME, 'Data/rRNA/only_mod.bed')
df_mod = pd.read_csv(mod_file, names=['sample', 'start', 'stop', 'mod'], sep='\t')
df_mod_Am = df_mod[(df_mod['mod']=='Am') & (df_mod['sample']=='NR_003286_RNA18SN5')]

ref_file = os.path.join(HOME, 'Data/transcriptomes/rRNA_18S_28S.fasta')
ref = {}
print('Parsing reference...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### WT ###
wt_bam_file = os.path.join(HOME, 'inference/rRNA/HCT116_wt_rRNA_sorted.bam')
wt_bam = pysam.AlignmentFile(wt_bam_file, 'rb')
# wt_fast5_dir = '/home/adrian/Data/rRNA/HCT116_wt_rRNA/fast5_pass'
wt_fast5_dir = '/prj/nanopore_direct_rnaseq/20210415_WFT_2/HCT116_wt_rRNA/20210415_1220_X2_AEV256_0b8f2ff4/fast5_pass'
wt_f5_paths = glob(os.path.join(wt_fast5_dir, '*.fast5'), recursive=True)
wt_index_read_ids = {}
print('Parsing WT fast5 files...', flush=True)
for f5_filepath in wt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        wt_index_read_ids[read_id] = f5_filepath
print('{} WT reads collected'.format(len(wt_index_read_ids)), flush=True)

### IVT ###
ivt_bam_file = os.path.join(HOME, 'inference/rRNA/18S_IVT_rRNA_sorted.bam')
ivt_bam = pysam.AlignmentFile(ivt_bam_file, 'rb')
# ivt_fast5_dir = '/home/adrian/Data/rRNA/18S_IVT_rRNA/fast5_pass/'
ivt_fast5_dir = '/prj/nanopore_direct_rnaseq/20210520_rRNA/18S_IVT_rRNA/20210520_1059_X1_AGT888_6ebcc4d8/fast5_pass'
ivt_f5_paths = glob(os.path.join(ivt_fast5_dir, '*.fast5'), recursive=True)
ivt_index_read_ids = {}
print('Parsing IVT fast5 files...', flush=True)
for f5_filepath in ivt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        ivt_index_read_ids[read_id] = f5_filepath
print('{} IVT reads collected'.format(len(ivt_index_read_ids)), flush=True)

### load model, device ###
model_path = os.path.join(HOME, 'pytorch_models/rRNA/rRNA-epoch29.torch')
torchdict = torch.load(model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(model_path, fixed_config)

### search by GLORI sites ###
MIN_COVERAGE = 0
outfile = os.path.join(HOME, 'Data/MAFIA/rRNA_outlier_ratios_thresh{:.2f}.tsv'.format(PERC_THRESH))
df_outlier = pd.DataFrame()
counts = 0
for ind, row in df_mod_Am.iterrows():
    ind = 0
    row = df_mod_Am.iloc[ind]
    print('\nSite {}'.format(ind), flush=True)
    sample = row['sample']
    site = row['start']
    mod = row['mod']
    ref_motif = ref[sample][site-2:site+3]

    # tic = time()
    wt_site_motif_features = collect_features_from_aligned_site(fixed_model, fixed_device, fixed_config, wt_bam, wt_index_read_ids, sample, site, MIN_COVERAGE)
    # elapsed = time() - tic
    # print('v2 time elapsed: {:.1f} secs'.format(elapsed))

    # print('Collecting IVT features...', flush=True)
    ivt_site_motif_features = collect_features_from_aligned_site_v2(fixed_model, fixed_device, fixed_config, ivt_bam, ivt_index_read_ids, chr, site, MIN_COVERAGE, enforce_motif=ref_motif)

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
