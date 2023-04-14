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
from unsupervised import train_cluster

parser = argparse.ArgumentParser()
parser.add_argument('--unm_bam_file')
parser.add_argument('--unm_fast5_dir')
parser.add_argument('--mod_bam_file')
parser.add_argument('--mod_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--mod_file')
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
mod_file = args.mod_file
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
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)
print('Parsed reference with contigs:\n', list(ref.keys()), flush=True)

df_mod = pd.read_excel(mod_file, skiprows=3)
print('Imported mod file with {} sites'.format(len(df_mod)), flush=True)

### unm ###
unm_bam = pysam.AlignmentFile(unm_bam_file, 'rb')
unm_f5_paths = glob(os.path.join(unm_fast5_dir, '*.fast5'), recursive=True)
unm_index_read_ids = {}
for f5_filepath in unm_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        unm_index_read_ids[read_id] = f5_filepath
print('{} IVT reads indexed'.format(len(unm_index_read_ids)), flush=True)

### mod ###
mod_bam = pysam.AlignmentFile(mod_bam_file, 'rb')
mod_f5_paths = glob(os.path.join(mod_fast5_dir, '*.fast5'), recursive=True)
mod_index_read_ids = {}
for f5_filepath in mod_f5_paths:
    f5 = get_fast5_file(f5_filepath, mode="r")
    for read_id in f5.get_read_ids():
        mod_index_read_ids[read_id] = f5_filepath
print('{} WT reads indexed'.format(len(mod_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(backbone_model_path, fixed_config, extraction_layer)

### loop through BID-seq sites ###
for ind, row in df_mod.iterrows():
    chr = row['chr'].lstrip('chr')
    pos = row['pos']
    name = row['name']
    strand = row['strand']
    bidseq_motif = row['Motif_1']
    bidseq_ratio = row['Frac_Ave %']
    if strand=='+':
        ref_motif = ref[chr][pos-2:pos+3]
    elif strand=='-':
        pos -= 1   ### why???
        ref_motif = ref[chr][pos-2:pos+3]
        ref_motif = str(Seq(ref_motif).reverse_complement())
    site_name = 'chr{} pos{} {}\nBID-seq ratio {}%'.format(chr, pos, name, round(bidseq_ratio))
    print('\n=========================================================', flush=True)
    print(site_name, flush=True)
    if ref_motif!=bidseq_motif:
        print('Aligned motif {} =/= BID-seq motif {}!'.format(ref_motif, bidseq_motif), flush=True)
        continue

    unm_motif_features = collect_features_from_aligned_site_v2(fixed_model, fixed_device, fixed_config, extraction_layer, unm_bam, unm_index_read_ids, chr, pos, min_coverage, max_num_reads, enforce_motif=ref_motif)
    print('{} IVT feature vectors collected'.format(len(unm_motif_features)), flush=True)
    mod_motif_features = collect_features_from_aligned_site_v2(fixed_model, fixed_device, fixed_config, extraction_layer, mod_bam, mod_index_read_ids, chr, pos, min_coverage, max_num_reads)
    print('{} WT feature vectors collected'.format(len(mod_motif_features)), flush=True)

    ### train classifier ###
    if (len(unm_motif_features)>=min_coverage) and (len(mod_motif_features)>=min_coverage):
        train_cluster(unm_motif_features, mod_motif_features, site_name, scaler, debug_img_dir=clustering_model_dir)
