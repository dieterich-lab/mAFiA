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
from extract_features import get_features_from_collection_of_signals, collect_site_features
from cluster_features import get_outlier_ratio_from_features, get_outlier_ratio_from_features_v2, get_outlier_ratio_from_features_v3
from time import time
import random
random.seed(10)
from random import sample
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--wt_bam_file', default=os.path.join(HOME, 'inference/rRNA/HCT116_wt_rRNA_sorted.bam'))
parser.add_argument('--wt_fast5_dir', default=os.path.join(HOME, 'Data/rRNA/HCT116_wt_rRNA/fast5_pass'))
parser.add_argument('--ivt_bam_file', default=os.path.join(HOME, 'inference/rRNA/18S_IVT_rRNA_sorted.bam'))
parser.add_argument('--ivt_fast5_dir', default=os.path.join(HOME, 'Data/rRNA/18S_IVT_rRNA/fast5_pass'))
parser.add_argument('--ref_file', default=os.path.join(HOME, 'Data/transcriptomes/rRNA_18S_28S.fasta'))
parser.add_argument('--mod_file', default=os.path.join(HOME, 'Data/rRNA/only_mod.bed'))
parser.add_argument('--rRNA_species', default='NR_003286_RNA18SN5')
parser.add_argument('--max_num_reads', default=100)
parser.add_argument('--min_coverage', default=0)
parser.add_argument('--mod_type', nargs='*', default=None, help='mod type')
parser.add_argument('--model_path', default=os.path.join(HOME, 'pytorch_models/rRNA/rRNA-epoch29.torch'))
parser.add_argument('--cluster_thresh', default=1.0, help='sigma threshold for clustering')
parser.add_argument('--outfile', default=os.path.join(HOME, 'inference/rRNA/mafia_outliers.tsv'))
args = parser.parse_args()
wt_bam_file = args.wt_bam_file
wt_fast5_dir = args.wt_fast5_dir
ivt_bam_file = args.ivt_bam_file
ivt_fast5_dir = args.ivt_fast5_dir
ref_file = args.ref_file
mod_file = args.mod_file
rRNA_species = args.rRNA_species
max_num_reads = int(args.max_num_reads)
min_coverage = int(args.min_coverage)
mod_type = args.mod_type
model_path = args.model_path
cluster_thresh = float(args.cluster_thresh)
outfile = os.path.join(HOME, 'inference/rRNA/{}_outlier_ratios_{}_sigma{:.2f}.tsv'.format(rRNA_species, mod_type, cluster_thresh))

df_mod = pd.read_csv(mod_file, names=['sample', 'start', 'stop', 'mod'], sep='\t')

if (mod_type is None) or (mod_type==[]):
    print('All mod types')
    mod_type = 'all'
    df_mod_sel = df_mod[df_mod['sample'] == rRNA_species]
else:
    print('Mod types: {}'.format(mod_type))
    df_mod_sel = df_mod[(df_mod['mod'].isin(mod_type)) & (df_mod['sample'] == rRNA_species)]

ref = {}
print('Parsing reference...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### WT ###
wt_bam = pysam.AlignmentFile(wt_bam_file, 'rb')
wt_f5_paths = glob(os.path.join(wt_fast5_dir, '*.fast5'), recursive=True)
wt_index_read_ids = {}
print('Parsing WT fast5 files...', flush=True)
for f5_filepath in wt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        wt_index_read_ids[read_id] = f5_filepath
print('{} WT reads collected'.format(len(wt_index_read_ids)), flush=True)

### IVT ###
ivt_bam = pysam.AlignmentFile(ivt_bam_file, 'rb')
ivt_f5_paths = glob(os.path.join(ivt_fast5_dir, '*.fast5'), recursive=True)
ivt_index_read_ids = {}
print('Parsing IVT fast5 files...', flush=True)
for f5_filepath in ivt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        ivt_index_read_ids[read_id] = f5_filepath
print('{} IVT reads collected'.format(len(ivt_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(model_path, fixed_config)

### extract features ###
print('Now extracting features from WT...')
wt_index_read_ids_sample = {id: wt_index_read_ids[id] for id in sample(list(wt_index_read_ids.keys()), min(len(wt_index_read_ids.keys()), max_num_reads))}
wt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, wt_index_read_ids_sample)

print('Now extracting features from IVT...')
ivt_index_read_ids_sample = {id: ivt_index_read_ids[id] for id in sample(list(ivt_index_read_ids.keys()), min(len(ivt_index_read_ids.keys()), max_num_reads))}
ivt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, ivt_index_read_ids_sample)

### search by GLORI sites ###
df_outlier = pd.DataFrame()
counts = 0
for ind, row in df_mod_sel.iterrows():
    # ind = 0
    # row = df_mod_sel.loc[ind]
    print('\nSite {}'.format(ind), flush=True)
    sample = row['sample']
    site = row['start']
    mod = row['mod']
    ref_motif = ref[sample][site-2:site+3]

    wt_site_motif_features = collect_site_features(wt_bam, sample, site, wt_predStr_features)
    # print('{} features from WT'.format(len(wt_site_motif_features)))
    ivt_site_motif_features = collect_site_features(ivt_bam, sample, site, ivt_predStr_features, enforce_motif=ref_motif)
    # print('{} features from IVT'.format(len(ivt_site_motif_features)))

    if (len(wt_site_motif_features)>min_coverage) and (len(ivt_site_motif_features)>min_coverage):
        print('=========================================================', flush=True)
        print('{}, pos{}, mod {}'.format(sample, site, mod), flush=True)
        print('Reference motif {}'.format(ref_motif), flush=True)
        print('{} feature vectors collected from WT'.format(len(wt_site_motif_features)), flush=True)
        print('{} feature vectors collected from IVT'.format(len(ivt_site_motif_features)), flush=True)
        print('Now clustering features...', flush=True)
        outlier_ratio = get_outlier_ratio_from_features(ivt_site_motif_features, wt_site_motif_features, ref_motif, cluster_thresh)
        # outlier_ratio = get_outlier_ratio_from_features_v2(ivt_site_motif_features, wt_site_motif_features, ref_motif, cluster_thresh)
        # outlier_ratio = get_outlier_ratio_from_features_v3(ivt_site_motif_features, wt_site_motif_features, ref_motif, cluster_thresh)
        print('Calculated outlier ratio {:.2f}'.format(outlier_ratio), flush=True)
        print('=========================================================', flush=True)
        if outlier_ratio!=-1:
            new_row = row.copy()
            new_row['motif'] = ref_motif
            new_row['ratio_outlier'] = np.round(outlier_ratio, 2)
            new_row['num_features_ivt'] = len(ivt_site_motif_features)
            new_row['num_features_wt'] = len(wt_site_motif_features)
            df_outlier = pd.concat([df_outlier, new_row.to_frame().T])
            counts += 1
            if counts%5==0:
                df_outlier.to_csv(outfile, sep='\t')