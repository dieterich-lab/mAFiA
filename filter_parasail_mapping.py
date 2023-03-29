import pandas as pd
import pysam
import numpy as np
import argparse
from glob import glob
from tqdm import tqdm

# fasta_file = '/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/achan/basecall/tRNA_IVT.fasta'
# ref_file = '/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta'
# csv_file = '/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/achan/mapping/parasail/tRNA_IVT.csv'
# rand_csv_files = glob('/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/achan/mapping/parasail/tRNA_IVT.csv.rand*')
# sam_file = '/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/achan/mapping/parasail/tRNA_IVT.sam'

parser = argparse.ArgumentParser()
parser.add_argument('--fasta_file')
parser.add_argument('--ref_file')
parser.add_argument('--csv_file')
parser.add_argument('--sam_file')

args = parser.parse_args()
fasta_file = args.fasta_file
ref_file = args.ref_file
csv_file = args.csv_file
rand_csv_files = glob(args.csv_file + '.rand*')
sam_file = args.sam_file

df_csv = pd.read_csv(csv_file, sep=',', names=['query', 'ref', 'query_len', 'ref_len', 'score', 'query_end', 'ref_end'])
dfs_csv_rand = [pd.read_csv(rand_csv_file, sep=',', names=['query', 'ref', 'query_len', 'ref_len', 'score', 'query_end', 'ref_end']) for rand_csv_file in rand_csv_files]

print('Parsing fasta file...')
query_names = []
with open(fasta_file, 'r') as f:
    for l in f.readlines():
        if l[0]=='>':
            query_names.append(l.lstrip('>').rstrip('\n'))

print('Parsing reference file...')
ref_names = []
with open(ref_file, 'r') as f:
    for l in f.readlines():
        if l[0]=='>':
            ref_names.append(l.lstrip('>').rstrip('\n'))

### calculate norm score ###
rand_scores = np.vstack([df_csv_rand['score'].values for df_csv_rand in dfs_csv_rand]).T
mu = np.mean(rand_scores, axis=1)
sigma = np.std(rand_scores, axis=1)
df_csv['norm_score'] = (df_csv['score'] - mu) / sigma

print('Matching query to reference...')
query_ref_matches = {}
for query_ind in tqdm(df_csv['query'].unique()):
    sub_df = df_csv[df_csv['query']==query_ind]
    ref_ind = sub_df.iloc[sub_df['norm_score'].argmax()]['ref']
    query_ref_matches[query_names[query_ind]] = ref_names[ref_ind]

print('Writing out filtered reads...')
with pysam.AlignmentFile(sam_file, 'r') as in_sam:
    with pysam.AlignmentFile(sam_file+'.filtered', 'w', template=in_sam) as out_sam:
        for read in in_sam:
            qname = read.query_name
            rname = read.reference_name
            if query_ref_matches[qname]==rname:
                out_sam.write(read)
