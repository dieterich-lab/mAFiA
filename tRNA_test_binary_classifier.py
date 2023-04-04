import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import argparse
from glob import glob
import pandas as pd
import torch
from models import objectview
import pysam
from Bio import SeqIO
from utils import index_fast5_files
from extract_features import load_model
from extract_features import get_features_from_collection_of_signals, collect_site_features
# from extract_features import collect_features_from_aligned_site_v2
from feature_classifiers import get_mod_ratio_with_binary_classifier
import random
random.seed(10)
from random import sample
from joblib import load

parser = argparse.ArgumentParser()
parser.add_argument('--test_bam_file')
parser.add_argument('--test_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--mod_file')
parser.add_argument('--max_num_reads', default=-1)
parser.add_argument('--min_coverage', default=0)
parser.add_argument('--backbone_model_path')
parser.add_argument('--extraction_layer')
parser.add_argument('--feature_width', default=0)
parser.add_argument('--classifier')
parser.add_argument('--classifier_model_dir')
parser.add_argument('--mod_prob_thresh', default=0)
parser.add_argument('--outfile')

args = parser.parse_args()
test_bam_file = args.test_bam_file
test_fast5_dir = args.test_fast5_dir
ref_file = args.ref_file
mod_file = args.mod_file
max_num_reads = int(args.max_num_reads)
min_coverage = int(args.min_coverage)
backbone_model_path = args.backbone_model_path
extraction_layer = args.extraction_layer
feature_width = int(args.feature_width)
classifier = args.classifier
classifier_model_dir = args.classifier_model_dir
mod_prob_thresh = float(args.mod_prob_thresh)
outfile = args.outfile

outdir = os.path.dirname(outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

df_mod = pd.read_csv(mod_file, names=['contig', 'start', 'stop', 'mod'], sep='\t')

ref = {}
print('Parsing genome...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### collect reads ###
test_bam = pysam.AlignmentFile(test_bam_file, 'rb')
test_f5_paths = glob(os.path.join(test_fast5_dir, '*.fast5'), recursive=True)
print('Parsing test fast5 files in {}...'.format(test_fast5_dir), flush=True)
test_index_read_ids = index_fast5_files(test_f5_paths)
print('{} test reads indexed'.format(len(test_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(backbone_model_path, fixed_config, extraction_layer)

### extract features ###
if max_num_reads>0:
    print('Now extracting features from WT...')
    test_index_read_ids_sample = {id: test_index_read_ids[id] for id in sample(list(test_index_read_ids.keys()), min(len(test_index_read_ids.keys()), max_num_reads))}
    test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids_sample, extraction_layer, feature_width)
else:
    print('Now extracting features from WT...')
    test_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, test_index_read_ids, extraction_layer, feature_width)

### loop through sites ###
df_out = pd.DataFrame()
counts = 0
for ind, mod_site in df_mod.iterrows():
    print('\nSite {}'.format(ind), flush=True)
    contig = mod_site['contig']
    start = int(mod_site['start'])
    mod = mod_site['mod']
    ref_motif = ref[contig][start - 2:start + 3]
    print('\n=========================================================', flush=True)
    print('{}, pos{}'.format(contig, start), flush=True)
    print('Reference motif {}'.format(ref_motif), flush=True)

    classifier_model_path = os.path.join(classifier_model_dir, '{}_{}_{}_{}.joblib'.format(classifier, contig, start, mod))
    if os.path.exists(classifier_model_path):
        classifier_model = load(classifier_model_path)
    else:
        print('Classifier model not found!\n', flush=True)
        print('=========================================================', flush=True)
        continue

    test_site_motif_features = collect_site_features(test_bam, contig, start, test_predStr_features)
    print('{} feature vectors collected'.format(len(test_site_motif_features)), flush=True)

    if len(test_site_motif_features)>min_coverage:
        clf = load(os.path.join(classifier_model_dir, '{}_{}_{}_{}.joblib'.format(classifier, contig, start, mod)))
        mod_ratio = get_mod_ratio_with_binary_classifier(test_site_motif_features, classifier_model, mod_prob_thresh)
        print('Predicted mod ratio {:.2f}'.format(mod_ratio), flush=True)
        print('=========================================================', flush=True)
        new_row = mod_site.copy()
        new_row['motif'] = ref_motif
        new_row['mod_ratio'] = round(mod_ratio, 3)
        new_row['num_test_features'] = len(test_site_motif_features)
        df_out = pd.concat([df_out, new_row.to_frame().T])
        counts += 1
        if counts % 5 == 0:
            df_out.to_csv(outfile, sep='\t')
df_out.to_csv(outfile, sep='\t')
print('Total {} sites'.format(counts))