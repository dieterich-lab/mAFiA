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
from extract_features import load_model, get_features_from_collection_of_signals, collect_all_motif_features
from feature_classifiers import get_mod_ratio_with_binary_classifier
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
parser.add_argument('--outfile')
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

### build dict of read and ref ###
dict_read_ref = {}
for read in test_bam.fetch():
    dict_read_ref[read.query_name] = read.reference_name

### load backbon model, device ###
torchdict = torch.load(args.backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(args.backbone_model_path, fixed_config, args.extraction_layer)

### load classifiers ###
classifier_model_paths = glob(os.path.join(args.classifier_model_dir, '*.joblib'))
classifier_models = {os.path.basename(this_path).rstrip('.joblib').split('_')[-1] : load(this_path) for this_path in classifier_model_paths}
classifier_motifs = list(classifier_models.keys())
print('Target motifs: {}'.format(', '.join(classifier_motifs)), flush=True)

### extract features ###
print('Now extracting features from test...', flush=True)
if args.max_num_reads>0:
    test_index_read_ids_sample = {id: test_index_read_ids[id] for id in sample(list(test_index_read_ids.keys()), min(len(test_index_read_ids.keys()), args.max_num_reads))}
    test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids_sample, args.extraction_layer, args.feature_width)
else:
    test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids, args.extraction_layer, args.feature_width)

### build dict of motif index and block size ###
index_motif_size_center = []
for k in ref.keys():
    if k.lstrip('block').split('_')[0] == '1':
        block_index = k.split('_')[1]
        block_seq = ref[k]
        block_size = len(block_seq)
        block_center = block_size // 2
        motif = block_seq[block_center-2:block_center+3]
        index_motif_size_center.append((block_index, motif, block_size, block_center))

### generate test results ###
outdir = os.path.dirname(args.outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

df_out = pd.DataFrame([])
for motif_ind, motif, block_size, block_center in index_motif_size_center:
    if motif not in classifier_motifs:
        print('Classifier for {} not available. Skipping...'.format(motif), flush=True)
        continue
    print('Now collecting features for motif {} from test reads...'.format(motif), flush=True)
    if args.enforce_motif:
        readId_qPos_motif_features = collect_all_motif_features(motif_ind, ref, test_bam, test_predStr_features, block_size=block_size, block_center=block_center, enforce_motif=motif)
    else:
        readId_qPos_motif_features = collect_all_motif_features(motif_ind, ref, test_bam, test_predStr_features, block_size=block_size, block_center=block_center)
    print('{} feature vectors collected'.format(len(readId_qPos_motif_features)), flush=True)

    if len(readId_qPos_motif_features)>args.min_coverage:
        _, readId_qPos_predMotif_modProbs = get_mod_ratio_with_binary_classifier(readId_qPos_motif_features, classifier_models[motif], output_mod_probs=True)
        df_out = pd.concat([
            df_out,
            pd.DataFrame([(x[0], dict_read_ref[x[0]], x[1], motif, x[2], round(x[3], 3)) for x in readId_qPos_predMotif_modProbs],
                         columns=['read_id', 'ref_name', 'q_pos', 'ref_motif', 'pred_motif', 'mod_prob'])
        ])
    df_out.to_csv(args.outfile, sep='\t', index=False)
print('Total number of bases tested {}'.format(len(df_out)), flush=True)