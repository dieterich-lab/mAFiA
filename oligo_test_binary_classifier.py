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
from feature_extractors import load_model, get_features_from_collection_of_signals, get_single_motif_nucleotides
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
parser.add_argument('--enforce_ref_5mer', action='store_true')
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
test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids, args.max_num_reads, args.extraction_layer, args.feature_width)

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
    print('Now collecting nucleotides for motif {}'.format(motif), flush=True)
    this_motif_nts = get_single_motif_nucleotides(motif_ind, ref, test_bam, test_predStr_features, block_size=block_size, block_center=block_center, enforce_ref_5mer=args.enforce_ref_5mer)
    print('{} NTs collected'.format(len(this_motif_nts)), flush=True)

    _ = get_mod_ratio_with_binary_classifier(this_motif_nts, classifier_models[motif])

    df_out = pd.concat([
        df_out,
        pd.DataFrame([(nt.read_id, dict_read_ref[nt.read_id], nt.read_pos, nt.ref_pos, motif, nt.pred_5mer, round(nt.mod_prob, 3)) for nt in this_motif_nts],
                     columns=['read_id', 'contig', 'q_pos', 't_pos', 'ref_motif', 'pred_motif', 'mod_prob'])
    ])
    df_out.to_csv(args.outfile, sep='\t', index=False)
print('Total number of nucleotides tested {}'.format(len(df_out)), flush=True)