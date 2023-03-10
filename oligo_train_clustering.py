import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import argparse
from glob import glob
import pandas as pd
import numpy as np
import torch
from models import objectview
import pysam
from Bio import SeqIO
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import load_model
from extract_features import get_features_from_collection_of_signals, collect_all_motif_features
# from feature_classifiers import train_binary_classifier
from unsupervised import train_cluster
import random
random.seed(10)
from random import sample
from joblib import dump

parser = argparse.ArgumentParser()
parser.add_argument('--unm_bam_file')
parser.add_argument('--unm_fast5_dir')
parser.add_argument('--mod_bam_file')
parser.add_argument('--mod_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--max_num_reads', default=-1)
parser.add_argument('--min_coverage', default=0)
parser.add_argument('--backbone_model_path')
parser.add_argument('--extraction_layer', default='convlayers.conv21')
parser.add_argument('--feature_width', default=0)
parser.add_argument('--scaler')
parser.add_argument('--clustering_model_dir')

args = parser.parse_args()
unm_bam_file = args.unm_bam_file
unm_fast5_dir = args.unm_fast5_dir
mod_bam_file = args.mod_bam_file
mod_fast5_dir = args.mod_fast5_dir
ref_file = args.ref_file
max_num_reads = int(args.max_num_reads)
min_coverage = int(args.min_coverage)
backbone_model_path = args.backbone_model_path
extraction_layer = args.extraction_layer
feature_width = int(args.feature_width)
scaler = args.scaler
clustering_model_dir = args.clustering_model_dir

if clustering_model_dir is not None:
    os.makedirs(clustering_model_dir, exist_ok=True)

ref = {}
print('Parsing reference...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### unm ###
unm_bam = pysam.AlignmentFile(unm_bam_file, 'rb')
unm_f5_paths = glob(os.path.join(unm_fast5_dir, '*.fast5'), recursive=True)
unm_index_read_ids = {}
print('Parsing unm fast5 files...', flush=True)
for f5_filepath in unm_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        unm_index_read_ids[read_id] = f5_filepath
print('{} unm reads indexed'.format(len(unm_index_read_ids)), flush=True)

### mod ###
mod_bam = pysam.AlignmentFile(mod_bam_file, 'rb')
mod_f5_paths = glob(os.path.join(mod_fast5_dir, '*.fast5'), recursive=True)
mod_index_read_ids = {}
print('Parsing mod fast5 files...', flush=True)
for f5_filepath in mod_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        mod_index_read_ids[read_id] = f5_filepath
print('{} mod reads indexed'.format(len(mod_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(backbone_model_path, fixed_config, extraction_layer)

### extract features ###
if max_num_reads>0:
    print('Now extracting features from unm...', flush=True)
    unm_index_read_ids_sample = {id: unm_index_read_ids[id] for id in sample(list(unm_index_read_ids.keys()), min(len(unm_index_read_ids.keys()), max_num_reads))}
    unm_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, unm_index_read_ids_sample, extraction_layer, feature_width)
    print('Now extracting features from mod...', flush=True)
    mod_index_read_ids_sample = {id: mod_index_read_ids[id] for id in sample(list(mod_index_read_ids.keys()), min(len(mod_index_read_ids.keys()), max_num_reads))}
    mod_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, mod_index_read_ids_sample, extraction_layer, feature_width)
else:
    print('Now extracting features from unm...', flush=True)
    unm_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, unm_index_read_ids, extraction_layer, feature_width)
    print('Now extracting features from mod...', flush=True)
    mod_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, mod_index_read_ids, extraction_layer, feature_width)

### collect motif features ###
# motif_ind = '0'
# motif_name = 'GGACA'

motif_indices_names = [
    ('0', 'GGACA'),
    ('1', 'GGACC'),
    ('2', 'AGACT')
]

for motif_ind, motif_name in motif_indices_names:
    print('Now collecting features for motif {} from unm reads...'.format(motif_name), flush=True)
    unm_motif_features = collect_all_motif_features(motif_ind, ref, unm_bam, unm_predStr_features, enforce_motif=True)
    print('{} feature vectors collected'.format(len(unm_motif_features)), flush=True)
    print('Now collecting features for motif {} from mod reads...'.format(motif_name), flush=True)
    mod_motif_features = collect_all_motif_features(motif_ind, ref, mod_bam, mod_predStr_features)
    print('{} feature vectors collected'.format(len(mod_motif_features)), flush=True)

    ### train classifier ###
    train_cluster(unm_motif_features, mod_motif_features, motif_name, scaler, debug_img_dir=os.path.join(clustering_model_dir, 'clustering'))

    # dump(classifier_model, os.path.join(classifier_model_dir, '{}_{}.joblib'.format(classifier, motif_name)))
    # print('AUC {:.2f}'.format(auc_score), flush=True)
