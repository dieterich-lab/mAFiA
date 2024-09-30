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
out_dir = os.path.join(res_dir, 'transcript_stoichiometry')
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


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}


def get_mod_probs_per_gene(in_bam, in_df):
    gene_mod_probs = {}
    for _, this_row in tqdm(in_df.iterrows()):
        chrom, chromStart, chromEnd, strand, gene = this_row[['chrom', 'chromStart', 'chromEnd', 'strand', 'gene']]
        gene_mod_probs[gene] = {mod: [] for mod in mods}
        flag_required = 0 if strand == '+' else 16
        for read in in_bam.fetch(chrom, chromStart, chromEnd):
            if read.flag == flag_required:
                for mod in mods:
                    gene_mod_probs[gene][mod].extend([pair[1] for pair in read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])])
    return gene_mod_probs


def get_num_nts_log2fc_pval(in_mod_probs_0, in_mod_probs_1):
    num_nts_0 = len(in_mod_probs_0)
    num_nts_1 = len(in_mod_probs_1)
    if (num_nts_0 > 0) and (num_nts_1 > 0):
        stoichio_0 = np.mean(np.array(in_mod_probs_0) >= 128)
        stoichio_1 = np.mean(np.array(in_mod_probs_1) >= 128)
        log2fc_stoichio = np.log2(stoichio_1 / stoichio_0)
        ks_stat, pval = kstest(in_mod_probs_0, in_mod_probs_1)
        return num_nts_0, num_nts_1, stoichio_0, stoichio_1, log2fc_stoichio, pval
    else:
        return None

with pysam.AlignmentFile(os.path.join(res_dir, ds, conditions[0], 'chrALL.mAFiA.reads.bam'), 'rb') as bam0:
    mod_probs_0 = get_mod_probs_per_gene(bam0, df_gene)

with pysam.AlignmentFile(os.path.join(res_dir, ds, conditions[1], 'chrALL.mAFiA.reads.bam'), 'rb') as bam1:
    mod_probs_1 = get_mod_probs_per_gene(bam1, df_gene)


gene_mod_log2fc = {}
for this_gene in tqdm(df_gene['gene']):
    gene_mod_log2fc[this_gene] = {}
    for this_mod in mods:
        gene_mod_log2fc[this_gene][this_mod] = get_num_nts_log2fc_pval(mod_probs_0[this_gene][this_mod],
                                                                       mod_probs_1[this_gene][this_mod])
gene_mod_log2fc = {k: v for k, v in gene_mod_log2fc.items() if any([vv is not None for vv in v.values()])}

df_out = pd.DataFrame([[k] + [kk] + list(vv) for k, v in gene_mod_log2fc.items() for kk, vv in v.items() if vv is not None],
                      columns=['gene',
                               'mod',
                               f'num_nts_{conditions[0]}',
                               f'num_nts_{conditions[1]}',
                               f'stoichiometry_{conditions[0]}',
                               f'stoichiometry_{conditions[1]}',
                               'log2fc',
                               'pval'
                               ]
                      )

df_out.to_csv(os.path.join(out_dir, f'{ds}_{conditions[1]}_vs_{conditions[0]}.tsv'),
              sep='\t', index=False, float_format='%.3f')
