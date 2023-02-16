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
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import load_model
from extract_features import get_features_from_collection_of_signals, collect_site_features
from cluster_features import train_svm_ivt_wt
import random
random.seed(10)
from random import sample
from joblib import dump

parser = argparse.ArgumentParser()
parser.add_argument('--wt_bam_file', default=os.path.join(HOME, 'inference/rRNA/HCT116_wt_rRNA_sorted.bam'))
parser.add_argument('--wt_fast5_dir', default=os.path.join(HOME, 'Data/rRNA/HCT116_wt_rRNA/fast5_pass'))
parser.add_argument('--ivt_bam_file', default=os.path.join(HOME, 'inference/rRNA/18S_IVT_rRNA_sorted.bam'))
parser.add_argument('--ivt_fast5_dir', default=os.path.join(HOME, 'Data/rRNA/18S_IVT_rRNA/fast5_pass'))
parser.add_argument('--ref_file', default=os.path.join(HOME, 'Data/transcriptomes/rRNA_18S_28S.fasta'))
parser.add_argument('--mod_file', default=os.path.join(HOME, 'Data/rRNA/only_mod.bed'))
parser.add_argument('--rRNA_species', default='NR_003286_RNA18SN5')
parser.add_argument('--max_num_reads', default=-1)
parser.add_argument('--min_coverage', default=0)
parser.add_argument('--mod_type', nargs='*', default=None, help='mod type')
parser.add_argument('--model_path', default=os.path.join(HOME, 'pytorch_models/rRNA/rRNA-epoch29.torch'))
parser.add_argument('--outfile', default=None)
parser.add_argument('--svm_model_dir', default=None)

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
svm_model_dir = args.svm_model_dir
if svm_model_dir is not None:
    os.makedirs(svm_model_dir, exist_ok=True)
outfile = args.outfile

df_mod = pd.read_csv(mod_file, names=['sample', 'start', 'stop', 'mod'], sep='\t')

if (mod_type is None) or (mod_type==[]):
    print('All mod types')
    df_mod_sel = df_mod[df_mod['sample'] == rRNA_species]
    if outfile is None:
        outfile = os.path.join(HOME, 'inference/rRNA/{}_mod_ratios_[{}].tsv'.format(rRNA_species, 'all'))
else:
    print('Mod types: {}'.format(mod_type))
    df_mod_sel = df_mod[(df_mod['mod'].isin(mod_type)) & (df_mod['sample'] == rRNA_species)]
    if outfile is None:
        outfile = os.path.join(HOME, 'inference/rRNA/{}_mod_ratios_[{}].tsv'.format(rRNA_species, '_'.join(mod_type)))

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
print('{} WT reads indexed'.format(len(wt_index_read_ids)), flush=True)

### IVT ###
ivt_bam = pysam.AlignmentFile(ivt_bam_file, 'rb')
ivt_f5_paths = glob(os.path.join(ivt_fast5_dir, '*.fast5'), recursive=True)
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
fixed_model, fixed_device = load_model(model_path, fixed_config)

### extract features ###
if max_num_reads>0:
    print('Now extracting features from WT...')
    wt_index_read_ids_sample = {id: wt_index_read_ids[id] for id in sample(list(wt_index_read_ids.keys()), min(len(wt_index_read_ids.keys()), max_num_reads))}
    wt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, wt_index_read_ids_sample)
    print('Now extracting features from IVT...')
    ivt_index_read_ids_sample = {id: ivt_index_read_ids[id] for id in sample(list(ivt_index_read_ids.keys()), min(len(ivt_index_read_ids.keys()), max_num_reads))}
    ivt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, ivt_index_read_ids_sample)
else:
    print('Now extracting features from WT...')
    wt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, wt_index_read_ids)
    print('Now extracting features from IVT...')
    ivt_predStr_features = get_features_from_collection_of_signals(fixed_model, fixed_device, fixed_config, ivt_index_read_ids)

df_mod_ratio = pd.DataFrame()
counts = 0
for ind, row in df_mod_sel.iterrows():
    print('\nSite {}'.format(ind), flush=True)
    sample = row['sample']
    start = int(row['start'])
    mod = row['mod']
    ref_motif = ref[sample][start - 2:start + 3]
    wt_site_motif_features = collect_site_features(wt_bam, sample, start, wt_predStr_features)
    ivt_site_motif_features = collect_site_features(ivt_bam, sample, start, ivt_predStr_features, enforce_motif=ref_motif)

    if (len(wt_site_motif_features)>min_coverage) and (len(ivt_site_motif_features)>min_coverage):
        print('=========================================================', flush=True)
        print('{}, pos{}, mod {}'.format(sample, start, mod), flush=True)
        print('Reference motif {}'.format(ref_motif), flush=True)
        print('{} feature vectors collected from WT'.format(len(wt_site_motif_features)), flush=True)
        print('{} feature vectors collected from IVT'.format(len(ivt_site_motif_features)), flush=True)
        print('Now classifying with SVM...', flush=True)
        acc, svm_model = train_svm_ivt_wt(ivt_site_motif_features, wt_site_motif_features, ref_motif)
        print('Accuracy {:.2f}'.format(acc), flush=True)
        print('=========================================================', flush=True)
        new_row = row.copy()
        new_row['motif'] = ref_motif
        new_row['svm_accuracy'] = int(acc*100)
        new_row['num_features_ivt'] = len(ivt_site_motif_features)
        new_row['num_features_wt'] = len(wt_site_motif_features)
        df_mod_ratio = pd.concat([df_mod_ratio, new_row.to_frame().T])
        if svm_model_dir is not None:
            dump(svm_model, os.path.join(svm_model_dir, 'svm_{}_{}_{}.joblib'.format(row['sample'], row['start'], row['mod'])))
        counts += 1
        if counts%5==0:
            df_mod_ratio.to_csv(outfile, sep='\t')
df_mod_ratio.to_csv(outfile, sep='\t')