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
from extract_features import get_features_from_collection_of_signals, collect_site_features
from feature_classifiers import train_svm_ivt_wt, train_binary_classifier
import random
random.seed(10)
from random import sample
from joblib import dump

parser = argparse.ArgumentParser()
parser.add_argument('--wt_bam_file')
parser.add_argument('--wt_fast5_dir')
parser.add_argument('--ivt_bam_file')
parser.add_argument('--ivt_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--mod_file')
parser.add_argument('--max_num_reads', default=-1)
parser.add_argument('--min_coverage', default=0)
parser.add_argument('--mod_type', nargs='*', default=None, help='mod type')
parser.add_argument('--model_path')
parser.add_argument('--extraction_layer', default='convlayers.conv21')
parser.add_argument('--feature_width', default=0)
parser.add_argument('--scaler', default=None)
parser.add_argument('--classifier')
parser.add_argument('--classifier_model_dir')
parser.add_argument('--outfile')

args = parser.parse_args()
wt_bam_file = args.wt_bam_file
wt_fast5_dir = args.wt_fast5_dir
ivt_bam_file = args.ivt_bam_file
ivt_fast5_dir = args.ivt_fast5_dir
ref_file = args.ref_file
mod_file = args.mod_file
max_num_reads = int(args.max_num_reads)
min_coverage = int(args.min_coverage)
mod_type = args.mod_type
model_path = args.model_path
extraction_layer = args.extraction_layer
feature_width = int(args.feature_width)
scaler = args.scaler
classifier = args.classifier
classifier_model_dir = args.classifier_model_dir
if classifier_model_dir is not None:
    os.makedirs(classifier_model_dir, exist_ok=True)
outfile = args.outfile

df_mod = pd.read_csv(mod_file, names=['contig', 'start', 'stop', 'mod'], sep='\t')

if (mod_type is None) or (mod_type==[]):
    print('All mod types')
    df_mod_sel = df_mod
else:
    print('Mod types: {}'.format(mod_type))
    df_mod_sel = df_mod[df_mod['mod'].isin(mod_type)]

ref = {}
print('Parsing reference...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

### WT ###
wt_bam = pysam.AlignmentFile(wt_bam_file, 'rb')
wt_f5_paths = glob(os.path.join(wt_fast5_dir, '**/*.fast5'), recursive=True)
wt_index_read_ids = {}
print('Parsing WT fast5 files...', flush=True)
for f5_filepath in wt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        wt_index_read_ids[read_id] = f5_filepath
print('{} WT reads indexed'.format(len(wt_index_read_ids)), flush=True)

### IVT ###
ivt_bam = pysam.AlignmentFile(ivt_bam_file, 'rb')
ivt_f5_paths = glob(os.path.join(ivt_fast5_dir, '**/*.fast5'), recursive=True)
ivt_index_read_ids = {}
print('Parsing IVT fast5 files...', flush=True)
for f5_filepath in ivt_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        ivt_index_read_ids[read_id] = f5_filepath
print('{} IVT reads indexed'.format(len(ivt_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(model_path, fixed_config, extraction_layer)

### extract features ###
if max_num_reads>0:
    print('Now extracting features from WT...')
    wt_index_read_ids_sample = {id: wt_index_read_ids[id] for id in sample(list(wt_index_read_ids.keys()), min(len(wt_index_read_ids.keys()), max_num_reads))}
    wt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, wt_index_read_ids_sample, extraction_layer, feature_width)
    print('Now extracting features from IVT...')
    ivt_index_read_ids_sample = {id: ivt_index_read_ids[id] for id in sample(list(ivt_index_read_ids.keys()), min(len(ivt_index_read_ids.keys()), max_num_reads))}
    ivt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, ivt_index_read_ids_sample, extraction_layer, feature_width)
else:
    print('Now extracting features from WT...')
    wt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, wt_index_read_ids, extraction_layer, feature_width)
    print('Now extracting features from IVT...')
    ivt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, ivt_index_read_ids, extraction_layer, feature_width)

df_mod_ratio = pd.DataFrame()
counts = 0
for ind, mod_site in df_mod_sel.iterrows():
    print('\nSite {}'.format(ind), flush=True)
    contig = mod_site['contig']
    start = int(mod_site['start'])
    mod = mod_site['mod']
    ref_motif = ref[contig][start - 2:start + 3]
    wt_site_motif_features = collect_site_features(wt_bam, contig, start, wt_predStr_features)
    ivt_site_motif_features = collect_site_features(ivt_bam, contig, start, ivt_predStr_features, enforce_motif=ref_motif)

    if (len(wt_site_motif_features)>min_coverage) and (len(ivt_site_motif_features)>min_coverage):
        print('=========================================================', flush=True)
        print('{}, pos{}, mod {}'.format(contig, start, mod), flush=True)
        print('Reference motif {}'.format(ref_motif), flush=True)
        print('{} feature vectors collected from WT'.format(len(wt_site_motif_features)), flush=True)
        print('{} feature vectors collected from IVT'.format(len(ivt_site_motif_features)), flush=True)
        if classifier=='svm':
            print('Now classifying with SVM...', flush=True)
            auc_score, classifier_model, opt_thresh = train_svm_ivt_wt(ivt_site_motif_features, wt_site_motif_features, ref_motif, classifier, debug_img_dir=os.path.join(classifier_model_dir, 'auc'))
        elif classifier=='logistic_regression':
            print('Now classifying with logistic regression...', flush=True)
            auc_score, classifier_model, opt_thresh = train_binary_classifier(ivt_site_motif_features, wt_site_motif_features, classifier, scaler=scaler, debug_img_path=os.path.join(classifier_model_dir, 'auc', '{}_{}_{}_{}.png'.format(classifier, contig, start, mod)))
        else:
            print('Classifier unspecified!')
            break

        print('AUC {:.2f}'.format(auc_score), flush=True)
        print('=========================================================', flush=True)
        new_row = mod_site.copy()
        new_row['motif'] = ref_motif
        new_row['auc_score'] = np.round(auc_score, 2)
        new_row['opt_thresh'] = np.round(opt_thresh, 2)
        new_row['num_features_ivt'] = len(ivt_site_motif_features)
        new_row['num_features_wt'] = len(wt_site_motif_features)
        df_mod_ratio = pd.concat([df_mod_ratio, new_row.to_frame().T])
        if classifier_model_dir is not None:
            dump(classifier_model, os.path.join(classifier_model_dir, '{}_{}_{}_{}.joblib'.format(classifier, contig, start, mod)))
        counts += 1
        if counts%5==0:
            df_mod_ratio.to_csv(outfile, sep='\t')
df_mod_ratio.to_csv(outfile, sep='\t')