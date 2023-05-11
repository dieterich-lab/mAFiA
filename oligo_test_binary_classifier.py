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
from utils import index_fast5_files
from extract_features import load_model
from extract_features import get_features_from_collection_of_signals
import random
random.seed(10)
from random import sample
from joblib import load

parser = argparse.ArgumentParser()
parser.add_argument('--test_bam_file')
parser.add_argument('--test_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--max_num_reads', type=int, default=-1)
parser.add_argument('--min_coverage', type=int, default=0)
parser.add_argument('--enforce_motif', action='store_true')
parser.add_argument('--backbone_model_path')
parser.add_argument('--extraction_layer', default='convlayers.conv21')
parser.add_argument('--feature_width', type=int, default=0)
parser.add_argument('--classifier_model_dir')
args = parser.parse_args()

### load reference ###
ref = {}
print('Parsing reference...', flush=True)
for record in SeqIO.parse(args.ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### load bam and index fast5 reads ###
test_bam = pysam.AlignmentFile(args.test_bam_file, 'rb')
test_f5_paths = glob(os.path.join(args.test_fast5_dir, '*.fast5'), recursive=True)
print('Parsing test fast5 files from {}'.format(args.test_fast5_dir), flush=True)
test_index_read_ids = index_fast5_files(test_f5_paths, test_bam)
print('{} test reads indexed'.format(len(test_index_read_ids)), flush=True)

### load backbon model, device ###
torchdict = torch.load(args.backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(args.backbone_model_path, fixed_config, args.extraction_layer)

### load classifiers ###
classifier_model_paths = glob(os.path.join(args.classifier_model_dir, '*.joblib'))
classifier_models = {os.path.basename(this_path).rstrip('.joblib').split('_')[-1] : load(this_path) for this_path in classifier_model_paths}
target_motifs = list(classifier_models.keys())
print('Target motifs: {}'.format(', '.join(target_motifs)), flush=True)

### extract features ###
print('Now extracting features from test...', flush=True)
if args.max_num_reads>0:
    test_index_read_ids_sample = {id: test_index_read_ids[id] for id in sample(list(test_index_read_ids.keys()), min(len(test_index_read_ids.keys()), args.max_num_reads))}
    test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids_sample, args.extraction_layer, args.feature_width)
else:
    test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids, args.extraction_layer, args.feature_width)

oligo_ind = 0



def get_mod_prob_for_single_read(features, dict_classifiers):
    motif_prob = {}
    for motif, clf in dict_classifiers.items():
        motif_prob[motif] = clf.predict_proba(features)[:, 1]
    return motif_prob

for read_id, (predStr, read_features) in test_predStr_features.items():
    read_modProb = get_mod_prob_for_single_read(read_features, classifier_models)
    mask_a = np.array([int(x == 'A') for x in list(predStr)])

    motifs = np.array(list(read_modProb.keys()))
    mat_prob = np.vstack([prob * mask_a for prob in read_modProb.values()])

    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    plt.figure(figsize=(20, 6))
    plt.imshow(mat_prob, vmin=0.8)
    plt.yticks(range(len(motifs)), motifs, rotation=0)
    plt.xticks(np.arange(len(predStr)), predStr)