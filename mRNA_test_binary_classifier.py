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
from Bio.Seq import Seq
from ont_fast5_api.fast5_interface import get_fast5_file
from extract_features import load_model
from extract_features import collect_features_from_aligned_site_v2
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
parser.add_argument('--use_opt_thresh', default=False)
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
use_opt_thresh = bool(args.use_opt_thresh)
outfile = args.outfile

outdir = os.path.dirname(outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

df_mod = pd.read_csv(mod_file)

# ref_file = os.path.join(HOME, 'Data/genomes/GRCh38_96.fa')
ref = {}
print('Parsing genome...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

### collect reads ###
test_bam = pysam.AlignmentFile(test_bam_file, 'rb')
test_f5_paths = glob(os.path.join(test_fast5_dir, '*.fast5'), recursive=True)
test_index_read_ids = {}
print('Parsing test fast5 files...', flush=True)
for f5_filepath in test_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        test_index_read_ids[read_id] = f5_filepath
print('{} test reads indexed'.format(len(test_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(backbone_model_path, fixed_config, extraction_layer)

### loop through sites ###
target_motifs = ['GGACA', 'GGACC', 'AGACT']
classifier_models = {this_motif : load(os.path.join(classifier_model_dir, '{}_{}.joblib'.format(classifier, this_motif))) for this_motif in target_motifs}
# target_motif = 'GGACA'
# classifier_model = load(os.path.join(classifier_model_dir, '{}_{}.joblib'.format(classifier, target_motif)))
df_mod_ratio = pd.DataFrame()
counts = 0
for ind, row in df_mod.iterrows():
    # print('\nSite {}'.format(ind), flush=True)
    chr = row['Chr'].lstrip('chr')
    if (chr.isnumeric()==False) and (chr not in ['X', 'Y']):
        continue
    strand = row['Strand']
    start = row['Sites'] - 1   # 0-based
    glori_ratio = row['Ratio']
    ref_motif = ref[chr][start-2:start+3]
    if strand=='-':
        ref_motif = str(Seq(ref_motif).reverse_complement())

    if ref_motif not in target_motifs:
        continue

    test_site_motif_features = collect_features_from_aligned_site_v2(fixed_model, fixed_device, fixed_config, extraction_layer, test_bam, test_index_read_ids, chr, start, min_coverage)

    if len(test_site_motif_features)>min_coverage:
        print('\n=========================================================', flush=True)
        print('chr{}, pos{}'.format(chr, start), flush=True)
        print('Reference motif {}'.format(ref_motif), flush=True)
        print('{} feature vectors collected'.format(len(test_site_motif_features)), flush=True)
        mod_ratio = get_mod_ratio_with_binary_classifier(test_site_motif_features, classifier_models[ref_motif])
        print('Predicted mod ratio {:.2f} [GLORI {:.2f}]'.format(mod_ratio, glori_ratio), flush=True)
        print('=========================================================', flush=True)
        new_row = row.copy()
        new_row['motif'] = ref_motif
        new_row['mod_ratio'] = round(mod_ratio, 3)
        new_row['num_test_features'] = len(test_site_motif_features)
        df_mod_ratio = pd.concat([df_mod_ratio, new_row.to_frame().T])
        counts += 1
        if counts % 5 == 0:
            df_mod_ratio.to_csv(outfile, sep='\t')
df_mod_ratio.to_csv(outfile, sep='\t')