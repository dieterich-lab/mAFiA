import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import pysam
import numpy as np
import re
from scipy.stats import kstest
from tqdm import tqdm

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
out_dir = os.path.join(res_dir, 'transcript_abundance')
os.makedirs(out_dir, exist_ok=True)

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']

gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'frame',
    'description'
]
df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
# df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]


def get_coverage_per_gene(in_bam, in_df):
    gene_coverage = {}
    for _, this_row in tqdm(in_df.iterrows()):
        chrom, chromStart, chromEnd, strand, gene = this_row[['chrom', 'chromStart', 'chromEnd', 'strand', 'gene']]
        flag_required = 0 if strand == '+' else 16
        coverage = 0
        for read in in_bam.fetch(chrom, chromStart, chromEnd):
            if read.flag == flag_required:
                coverage += 1
        gene_coverage[gene] = coverage
    return gene_coverage


with pysam.AlignmentFile(os.path.join(res_dir, ds, conditions[0], 'chrALL.mAFiA.reads.bam'), 'rb') as bam0:
    coverage0 = get_coverage_per_gene(bam0, df_gene)

with pysam.AlignmentFile(os.path.join(res_dir, ds, conditions[1], 'chrALL.mAFiA.reads.bam'), 'rb') as bam1:
    coverage1 = get_coverage_per_gene(bam1, df_gene)

sum0 = np.sum(list(coverage0.values()))
sum1 = np.sum(list(coverage1.values()))

gene_coverage_log2fc = {}
for this_gene in df_gene['gene']:
    cov0 = coverage0[this_gene]
    cov1 = coverage1[this_gene]
    if (cov0 > 0) and (cov1 > 0):
        normCov0 = cov0 / sum0
        normCov1 = cov1 / sum1
        gene_coverage_log2fc[this_gene] = (cov0, cov1, normCov0, normCov1, np.log2(normCov1 / normCov0))

df_out = pd.DataFrame([[k] + list(v) for k, v in gene_coverage_log2fc.items()],
                      columns=['gene',
                               f'coverage_{conditions[0]}',
                               f'coverage_{conditions[1]}',
                               f'norm_coverage_{conditions[0]}',
                               f'norm_coverage_{conditions[1]}',
                               'log2fc'
                               ]
                      )
df_out[f'norm_coverage_{conditions[0]}'] = df_out[f'norm_coverage_{conditions[0]}'].map('{:.3E}'.format)
df_out[f'norm_coverage_{conditions[1]}'] = df_out[f'norm_coverage_{conditions[1]}'].map('{:.3E}'.format)
df_out[f'log2fc'] = df_out['log2fc'].map('{:.3f}'.format)

df_out.to_csv(os.path.join(out_dir, f'{ds}_{conditions[1]}_vs_{conditions[0]}.tsv'), sep='\t', index=False)
