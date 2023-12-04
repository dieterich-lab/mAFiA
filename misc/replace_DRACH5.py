import os
import pandas as pd

prj_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A'

MIN_COV = 50

ds = '100_WT_0_IVT'
# ds = '0_WT_100_IVT'
# ds = 'Mettl3-KO'
# chr = 14

for this_chr in list(range(1, 23)) + 'X':

    orig_file = os.path.join(prj_dir, f'DRACH/{ds}/chr{this_chr}/mAFiA.sites.bed')
    replace_file = os.path.join(prj_dir, f'DRACH5_retrain/{ds}/chr{this_chr}/mAFiA.sites.bed')
    out_dir = os.path.join(prj_dir, f'DRACH_replaced/{ds}/chr{this_chr}')

    os.makedirs(out_dir, exist_ok=True)
    df_orig = pd.read_csv(orig_file, sep='\t', dtype={'chrom': str})
    df_replace = pd.read_csv(replace_file, sep='\t', dtype={'chrom': str})
    replace_ref5mer = df_replace['ref5mer'].unique()

    df_new = pd.concat((
        df_orig[
            (~df_orig['ref5mer'].isin(replace_ref5mer))
            * (df_orig['coverage']>=MIN_COV)
        ],
        df_replace
    ), ignore_index=True)

    df_new.sort_values(by='chromStart', inplace=True)
    df_new.to_csv(os.path.join(out_dir, 'mAFiA.sites.bed'), sep='\t', index=False)