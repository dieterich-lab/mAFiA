import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import argparse
from glob import glob
import numpy as np
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# parser = argparse.ArgumentParser()
# parser.add_argument('--ref_file')
# parser.add_argument('--mod_file')
# parser.add_argument('--outfile')
# args = parser.parse_args()
# ref_file = args.ref_file
# mod_file = args.mod_file
# outfile = args.outfile

# COVERAGE_THRESH = 10
NUM_MOTIFS = 16
P_VAL_THRESH = 1.0E-10

ref_file = '/home/adrian/Data/genomes/GRCh38_96.fa'
mod_file = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv'
test_res_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_HEK293A_WT_multiple_NoNorm_allReads.tsv'
img_out = os.path.dirname(test_res_file)

outdir = os.path.dirname(test_res_file)
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

### GLORI mod file ###
df_mod = pd.read_csv(mod_file)

### parse reference ###
ref = {}
print('Parsing genome...', flush=True)
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
        ref[record.id] = str(record.seq)

all_ref_motifs = []
for ind, row in tqdm(df_mod.iterrows()):
    if row['Pvalue'] > P_VAL_THRESH:
        continue
    chr = row['Chr'].lstrip('chr')
    start = row['Sites'] - 1   # 0-based
    if (chr.isnumeric()==False) and (chr not in ['X', 'Y']):
        continue
    strand = row['Strand']
    glori_ratio = row['Ratio']
    ref_motif = ref[chr][start-2:start+3]
    if strand=='-':
        ref_motif = str(Seq(ref_motif).reverse_complement())
    all_ref_motifs.append(ref_motif)
# print('{} GLORI sites with p-val below {}'.format(len(all_ref_motifs), P_VAL_THRESH))

### test output ###
df_mafia = pd.read_csv(test_res_file, sep='\t', index_col=0)
# counts = len(df_mafia)
# print('Test result {} with {} sites'.format(os.path.basename(test_res_file), counts))

### histogram ###
glori_motif_counts = Counter(all_ref_motifs).most_common(NUM_MOTIFS)
percentage_glori = np.sum([v for k, v in glori_motif_counts]) / len(all_ref_motifs) * 100
glori_motif_keys = np.array([k for k, v in glori_motif_counts])

plt.figure(figsize=(NUM_MOTIFS, 6))
plt.bar(np.arange(len(glori_motif_counts)), [v for k, v in glori_motif_counts], label='GLORI m6A sites')
for cov_thresh in [10, 20, 50]:
    df_mafia_thresh = df_mafia[df_mafia['num_test_features']>=cov_thresh]
    mafia_motif_counts = Counter(df_mafia_thresh['motif']).most_common()
    plt.bar([np.where(glori_motif_keys==k)[0][0] for k, v in mafia_motif_counts], [v for k, v in mafia_motif_counts], label='ONT coverage$\geq${}'.format(cov_thresh))
plt.xticks(np.arange(len(glori_motif_counts)), glori_motif_keys, rotation=-90)
plt.legend(loc='upper right', fontsize=15)
plt.xlabel('Motif', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.title('First {} motifs\n{:.1f}% of {} GLORI sites (p-val$\leq${})'.format(NUM_MOTIFS, percentage_glori, len(all_ref_motifs), P_VAL_THRESH), fontsize=15)
plt.savefig(os.path.join(img_out, '{}_percentage_{}motifs_glori_pVal{}.png'.format(os.path.basename(test_res_file).replace('.tsv', ''), NUM_MOTIFS, P_VAL_THRESH)), bbox_inches='tight')
plt.close('all')