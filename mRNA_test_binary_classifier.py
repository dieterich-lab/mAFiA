import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import argparse
from glob import glob, iglob
from tqdm import tqdm
import pandas as pd
import torch
from models import objectview
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from utils import index_fast5_files
from extract_features import load_model
from extract_features import get_nucleotides_aligned_to_site
from feature_classifiers import get_mod_ratio_with_binary_classifier
from joblib import load

parser = argparse.ArgumentParser()
parser.add_argument('--test_bam_file')
parser.add_argument('--test_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--mod_file')
parser.add_argument('--max_num_reads', type=int, default=-1)
parser.add_argument('--min_coverage', type=int, default=0)
parser.add_argument('--enforce_ref_5mer', action='store_true')
parser.add_argument('--backbone_model_path')
parser.add_argument('--extraction_layer', default='convlayers.conv21')
parser.add_argument('--feature_width', type=int, default=0)
parser.add_argument('--classifier', default='logistic_regression')
parser.add_argument('--classifier_model_dir')
parser.add_argument('--mod_prob_thresh', type=float, default=0.5)
parser.add_argument('--outfile')
parser.add_argument('--output_mod_probs', action='store_true')
args = parser.parse_args()

outdir = os.path.dirname(args.outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

df_mod = pd.read_csv(args.mod_file)
df_mod = df_mod.rename(columns={'Unnamed: 0': 'index'})

ref = {}
print('Parsing genome reference {}...'.format(args.ref_file), flush=True)
for record in tqdm(SeqIO.parse(args.ref_file, 'fasta')):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

### collect reads ###
test_bam = pysam.AlignmentFile(args.test_bam_file, 'rb')
test_f5_paths = list(iglob(os.path.join(args.test_fast5_dir, '*.fast5'), recursive=True))
print('Parsing {} test fast5 files in {}'.format(len(test_f5_paths), args.test_fast5_dir), flush=True)
test_index_read_ids = index_fast5_files(test_f5_paths)
print('{} test reads indexed'.format(len(test_index_read_ids)), flush=True)

### load model, device ###
torchdict = torch.load(args.backbone_model_path, map_location="cpu")
origconfig = torchdict["config"]
fixed_config = objectview(origconfig)
fixed_model, fixed_device = load_model(args.backbone_model_path, fixed_config, args.extraction_layer)

### load classifiers ###
classifier_model_paths = glob(os.path.join(args.classifier_model_dir, '{}_*.joblib'.format(args.classifier)))
classifier_models = {os.path.basename(this_path).rstrip('.joblib').split('_')[-1] : load(this_path) for this_path in classifier_model_paths}
target_motifs = list(classifier_models.keys())
print('Target motifs: {}'.format(', '.join(target_motifs)), flush=True)

### check for existing output ###
restart = False
if os.path.exists(args.outfile):
    df_out = pd.read_csv(args.outfile, sep='\t')
    counts = len(df_out)
    if counts>0:
        print('Restarting from {} with {} counts'.format(args.outfile, counts), flush=True)
        last_row = df_out.tail(1)
        last_chr = last_row['Chr'].values[0].lstrip('chr')
        last_start = last_row['Sites'].values[0] - 1  # 0-based
        restart = True
else:
    df_out = pd.DataFrame()
    counts = 0
    print('Starting from scratch', flush=True)
    restart = False

### iterate over mod sites ###
for _, row in df_mod.iterrows():
    ind = row['index']
    chr = row['Chr'].lstrip('chr')
    start = row['Sites'] - 1   # 0-based

    if restart:
        if (chr==last_chr) and (start==last_start):
            print('Last row in file: chr{}, pos{}'.format(chr, start), flush=True)
            restart = False
        else:
            print('Skipping chr{}, pos{}'.format(chr, start), flush=True)
        continue

    if (chr.isnumeric()==False) and (chr not in ['X', 'Y']):
        continue
    strand = row['Strand']
    glori_ratio = row['Ratio']
    ref_motif = ref[chr][start-2:start+3]
    if strand=='-':
        ref_motif = str(Seq(ref_motif).reverse_complement())

    if ref_motif not in target_motifs:
        continue

    test_site_nucleotides = get_nucleotides_aligned_to_site(fixed_model, fixed_device, fixed_config, args.extraction_layer, test_bam, test_index_read_ids, chr, start, strand, args.min_coverage, args.max_num_reads, enforce_ref_5mer=args.enforce_ref_5mer, ref_5mer=ref_motif)

    if len(test_site_nucleotides)>args.min_coverage:
        print('=========================================================', flush=True)
        print('{}, chr{}, pos{}, strand{}'.format(ind, chr, start, strand), flush=True)
        print('Reference motif {}'.format(ref_motif), flush=True)
        print('{} feature vectors collected'.format(len(test_site_nucleotides)), flush=True)
        mod_ratio = get_mod_ratio_with_binary_classifier(test_site_nucleotides, classifier_models[ref_motif],
                                                         args.mod_prob_thresh)
        new_row = row.copy()
        new_row['ref_motif'] = ref_motif
        if args.output_mod_probs:
            for nt in test_site_nucleotides:
                new_row['pred_motif'] = nt.pred_5mer
                new_row['read_id'] = nt.read_id
                new_row['mod_prob'] = round(nt.mod_prob, 3)
                df_out = pd.concat([df_out, new_row.to_frame().T])
        else:
            new_row['mod_ratio'] = round(mod_ratio, 3)
            new_row['num_test_features'] = len(test_site_nucleotides)
            df_out = pd.concat([df_out, new_row.to_frame().T])
        print('Predicted mod ratio {:.2f} [GLORI {:.2f}]'.format(mod_ratio, glori_ratio), flush=True)
        print('=========================================================\n', flush=True)
        counts += 1
        df_out.to_csv(args.outfile, sep='\t', index=False)
print('Total {} sites written to {}'.format(counts, args.outfile), flush=True)