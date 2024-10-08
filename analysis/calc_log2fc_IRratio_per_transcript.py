import os
import pandas as pd
import numpy as np
from scipy.stats import kstest
from tqdm import tqdm

ds = 'TAC'
conditions = ['SHAM', 'TAC']

IRF_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/IRFinder'
cond_df = {
    this_cond: pd.read_csv(
        os.path.join(IRF_dir, f'{ds}_{this_cond}_merged.txt'), sep='\t', dtype={'Chr': str}
    ) for this_cond in conditions
}

df_merged = pd.merge(cond_df[conditions[0]], cond_df[conditions[1]],
                     on=['Chr', 'Start', 'End', 'Name', 'Null', 'Strand'],
                     suffixes=[f'_{conditions[0]}', f'_{conditions[1]}']
                     )
df_merged['gene'] = [this_name.split('/')[0] for this_name in df_merged['Name']]

unique_genes = list(df_merged['gene'].unique())
gene_log2fc_IRratio_pval = {}
for this_gene in tqdm(unique_genes):
    ratios_0 = df_merged[df_merged['gene'] == this_gene][f'IRratio_{conditions[0]}'].values
    ratios_1 = df_merged[df_merged['gene'] == this_gene][f'IRratio_{conditions[1]}'].values
    mean_0 = np.mean(ratios_0)
    mean_1 = np.mean(ratios_1)
    log2fc = np.log2(mean_1/mean_0)
    _, pval = kstest(ratios_0, ratios_1)
    gene_log2fc_IRratio_pval[this_gene] = (len(ratios_0), mean_0, mean_1, log2fc, pval)
df_gene_log2fc_IRratio_pval = pd.DataFrame([[k] + list(v) for k, v in gene_log2fc_IRratio_pval.items()
                                            if ~np.isnan(v[-2]) and ~np.isinf(v[-2])],
                                           columns=['gene',
                                                  'num_IR_junctions',
                                                  f'mean_IR_ratio_{conditions[0]}',
                                                  f'mean_IR_ratio_{conditions[0]}',
                                                  'log2fc',
                                                  'pval'
                                                  ])
df_gene_log2fc_IRratio_pval.sort_values('pval', inplace=True)
df_gene_log2fc_IRratio_pval.to_csv(os.path.join(IRF_dir, f'gene_log2fc_IRratio_pval_{ds}.tsv'), sep='\t', index=False)