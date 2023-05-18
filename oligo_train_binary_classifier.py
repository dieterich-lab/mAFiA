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
from utils import index_fast5_files
from extract_features import load_model
from extract_features import get_features_from_collection_of_signals, get_single_motif_nucleotides
from feature_classifiers import train_binary_classifier
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
parser.add_argument('--max_num_reads', type=int, default=-1)
parser.add_argument('--min_coverage', type=int, default=0)
parser.add_argument('--enforce_ref_5mer', action='store_true')
parser.add_argument('--backbone_model_path')
parser.add_argument('--extraction_layer', default='convlayers.conv21')
parser.add_argument('--feature_width', type=int, default=0)
parser.add_argument('--scaler', default=None)
parser.add_argument('--classifier', default='logistic_regression')
parser.add_argument('--classifier_model_dir')
args = parser.parse_args()

if args.enforce_ref_5mer:
    classifier_model_dir = args.classifier_model_dir + '_enforceMotif'
else:
    classifier_model_dir = args.classifier_model_dir
os.makedirs(classifier_model_dir, exist_ok=True)

ref = {}
print('Parsing reference...', flush=True)
for record in SeqIO.parse(args.ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### unm ###
unm_bam = pysam.AlignmentFile(args.unm_bam_file, 'rb')
unm_f5_paths = glob(os.path.join(args.unm_fast5_dir, '*.fast5'), recursive=True)
print('Parsing unm fast5 files from {}'.format(args.unm_fast5_dir), flush=True)
unm_index_read_ids = index_fast5_files(unm_f5_paths, unm_bam)
print('{} unm reads indexed'.format(len(unm_index_read_ids)), flush=True)

### mod ###
mod_bam = pysam.AlignmentFile(args.mod_bam_file, 'rb')
mod_f5_paths = glob(os.path.join(args.mod_fast5_dir, '*.fast5'), recursive=True)
print('Parsing mod fast5 files from {}'.format(args.mod_fast5_dir), flush=True)
mod_index_read_ids = index_fast5_files(mod_f5_paths, mod_bam)
print('{} mod reads indexed'.format(len(mod_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(args.backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(args.backbone_model_path, fixed_config, args.extraction_layer)

### extract features ###
print('Now extracting features from unm...', flush=True)
unm_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, unm_index_read_ids, args.max_num_reads, args.extraction_layer, args.feature_width)
print('Now extracting features from mod...', flush=True)
mod_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, mod_index_read_ids, args.max_num_reads, args.extraction_layer, args.feature_width)

### block center positions ###
block_index_motif_size_center = []
for k in ref.keys():
    if k.lstrip('block').split('_')[0] == '1':
        block_index = k.split('_')[1]
        block_seq = ref[k]
        block_size = len(block_seq)
        block_center = block_size // 2
        motif = block_seq[block_center-2:block_center+3]
        block_index_motif_size_center.append((block_index, motif, block_size, block_center))

for motif_ind, motif, block_size, block_center in block_index_motif_size_center:
    print('Now collecting nucleotides for motif {} from unm reads...'.format(motif), flush=True)
    unm_motif_nts = get_single_motif_nucleotides(motif_ind, ref, unm_bam, unm_predStr_features, block_size=block_size, block_center=block_center, enforce_ref_5mer=args.enforce_ref_5mer)
    print('{} NTs collected'.format(len(unm_motif_nts)), flush=True)
    print('Now collecting nucleotides for motif {} from mod reads...'.format(motif), flush=True)
    mod_motif_nts = get_single_motif_nucleotides(motif_ind, ref, mod_bam, mod_predStr_features, block_size=block_size, block_center=block_center, enforce_ref_5mer=args.enforce_ref_5mer)
    print('{} NTs collected'.format(len(mod_motif_nts)), flush=True)

    ### train classifier ###
    if min(len(unm_motif_nts), len(mod_motif_nts))>=args.min_coverage:
        auc_score, classifier_model, opt_thresh = train_binary_classifier(unm_motif_nts, mod_motif_nts, classifier=args.classifier, scaler=args.scaler, debug_img_path=os.path.join(classifier_model_dir, 'auc_{}.png'.format(motif)), fig_title=motif)
        dump(classifier_model, os.path.join(classifier_model_dir, '{}_{}.joblib'.format(args.classifier, motif)))
        print('AUC {:.2f}'.format(auc_score), flush=True)
