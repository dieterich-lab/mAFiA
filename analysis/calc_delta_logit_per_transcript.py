import os
import pandas as pd
import re
import numpy as np
from scipy.stats import kstest
from tqdm import tqdm

def get_logit(vec_modRatio):
    return np.log2(vec_modRatio / (100 - vec_modRatio))


########################################################################################################################
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

########################################################################################################################

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
out_dir = os.path.join(res_dir, 'transcript_logit')
os.makedirs(out_dir, exist_ok=True)

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']

# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

### merge conditions ###
cond_df = {}
for this_cond in conditions:
    res_file = os.path.join(res_dir, f'{ds}/{this_cond}/chrALL.mAFiA.sites.bed')
    cond_df[this_cond] = pd.read_csv(res_file, sep='\t', dtype={'chrom': str})

df_merged = pd.merge(*cond_df.values(),
                     on=['chrom',
                         'chromStart',
                         'chromEnd',
                         'name',
                         'score',
                         'strand',
                         'ref5mer'],
                     suffixes=[f"_{this_cond.rstrip('_merged')}" for this_cond in conditions])

### logit per transcript ###
gene_mod_sites_delta_pval = []
for this_gene in tqdm(df_gene['gene']):
    sel_gene_row = df_gene[df_gene['gene'] == this_gene].iloc[0]
    gene_df = df_merged[
        (df_merged.chrom == sel_gene_row.chrom)
        * (df_merged.chromStart >= sel_gene_row.chromStart)
        * (df_merged.chromEnd <= sel_gene_row.chromEnd)
        * (df_merged.strand == sel_gene_row.strand)
    ]
    for this_mod in ['m6A', 'psi']:
        mod_df = gene_df[gene_df.name == this_mod]

        logit1 = get_logit(mod_df[f"modRatio_{conditions[1].rstrip('_merged')}"].values)
        logit0 = get_logit(mod_df[f"modRatio_{conditions[0].rstrip('_merged')}"].values)

        vec_delta_logit = logit1 - logit0
        vec_delta_logit = vec_delta_logit[(~np.isinf(vec_delta_logit)) * (~np.isnan(vec_delta_logit))]

        if len(vec_delta_logit) == 0:
            continue

        xlim = [-2.5, 2.5]
        sigma = np.std(vec_delta_logit)
        vec_x = np.linspace(*xlim, 100)
        gaussian = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-vec_x**2 / (2 * sigma**2))

        ks_stat, pval = kstest(gaussian, vec_delta_logit)
        mean_delta = np.mean(vec_delta_logit)

        gene_mod_sites_delta_pval.append([this_gene, this_mod, len(vec_delta_logit), mean_delta, pval])


df_out = pd.DataFrame(gene_mod_sites_delta_pval, columns=['gene', 'mod', 'num_sites', 'delta_logit', 'pval'])
df_out.to_csv(os.path.join(out_dir, f'delta_logitS_{ds}.tsv'), sep='\t', index=False)