import os
from glob import glob
import pandas as pd
from collections import Counter
import numpy as np
from scipy.stats import kstest
from tqdm import tqdm


def get_df_log2fc(in_cond_df):
    common_genes = list(set(in_cond_df[conditions[0]]['gene'].unique()).intersection(set(in_cond_df[conditions[1]]['gene'].unique())))
    out_gene_polyA_log2fc_pval = {}
    for thisGene in tqdm(common_genes):
        polyA_len_0 = in_cond_df[conditions[0]][in_cond_df[conditions[0]]['gene'] == thisGene]['polyA_length']
        polyA_len_1 = in_cond_df[conditions[1]][in_cond_df[conditions[1]]['gene'] == thisGene]['polyA_length']
        mean_len_0 = np.mean(polyA_len_0)
        mean_len_1 = np.mean(polyA_len_1)
        log2fc = np.log2(mean_len_1 / mean_len_0)
        ks_stat, pval = kstest(polyA_len_0, polyA_len_1)
        out_gene_polyA_log2fc_pval[thisGene] = (len(polyA_len_0), len(polyA_len_1), mean_len_0, mean_len_1, log2fc, pval)

    out_df_gene_polyA_log2fc_pval = pd.DataFrame([[k] + list(v) for k, v in out_gene_polyA_log2fc_pval.items()],
                                             columns=['gene',
                                                      f'num_reads_{conditions[0]}',
                                                      f'num_reads_{conditions[1]}',
                                                      f'mean_polyA_len_{conditions[0]}',
                                                      f'mean_polyA_len_{conditions[1]}',
                                                      'log2fc',
                                                      'pval'
                                                      ])
    out_df_gene_polyA_log2fc_pval.sort_values('log2fc', inplace=True)
    return out_df_gene_polyA_log2fc_pval


polyA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA'

# ds = 'TAC'
# conditions = ['SHAM', 'TAC']
# ds = 'HFpEF'
# conditions = ['ctrl', 'HFpEF']
ds = 'Diet'
conditions = ['WT_CD', 'WT_WD']

log2fc_filename = os.path.join(polyA_dir, f'gene_polyA_log2fc_pval_{ds}.tsv')

cond_df = {}
for this_cond in conditions:
    all_filepaths = glob(os.path.join(polyA_dir, f'polyA_reads_annotated_{ds}_{this_cond}_merged.tsv'))
    all_dfs = []
    for this_filepath in all_filepaths:
        all_dfs.append(pd.read_csv(this_filepath, sep='\t'))
    cond_df[this_cond] = pd.concat(all_dfs).sort_values('gene')

if os.path.exists(log2fc_filename):
    df_gene_polyA_log2fc_pval = pd.read_csv(log2fc_filename, sep='\t')
else:
    df_gene_polyA_log2fc_pval = get_df_log2fc(cond_df)
    df_gene_polyA_log2fc_pval.to_csv(log2fc_filename, sep='\t', index=False, float_format='%.6f')
