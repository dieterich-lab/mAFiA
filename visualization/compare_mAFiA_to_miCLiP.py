import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt
import os
from collections import Counter

miCLiP_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

dict_threshold = {
    'confidence' : 80.0,
    'modRatio' : 50.0,
    'score_miCLiP' : 0.5
}

# ds = '40-34'
# cf = 'sham'

ds = '40-26'
cf = 'TAC'

mAFiA_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/mouse_heart/Dewenter_TAC_Backs_lab/{ds}/chrALL.mAFiA.sites.bed'
miCLiP_file = f'/home/adrian/Data/Heart_on_m6A/miCLIP/analysis_You_Zhou/nochr_predicted_m6A_{cf}.bed'
out_dir = '/home/adrian/Data/Heart_on_m6A/miCLIP/adrian'
out_df_file = os.path.join(out_dir, f'overlap_{ds}_cf_miCLiP_{cf}.tsv')
out_motif_file = os.path.join(out_dir, f'motifs_{ds}_cf_miCLiP_{cf}.tsv')

df_mAFiA = pd.read_csv(mAFiA_file, sep='\t', dtype={'chrom': str})
df_miCLiP = pd.read_csv(miCLiP_file, sep='\t', names=miCLiP_fields, dtype={'chrom': str})
df_miCLiP['name'] = 'm6A'

df_merged = pd.merge(df_miCLiP, df_mAFiA, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand'], suffixes=['_miCLiP', '_mAFiA'])
query = ' and '.join([f'{k} >= {repr(v)}' for k, v in dict_threshold.items()])
df_merged_sel = df_merged.query(query)

df_merged_sel = pd.concat([
    df_merged_sel[df_merged_sel['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_merged_sel[~df_merged_sel['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_merged_sel.drop(columns=['score_mAFiA'], inplace=True)
df_merged_sel.to_csv(out_df_file, sep='\t', index=False, float_format='%.3f')

with open(out_motif_file, 'w') as f_out:
    for k, v in Counter(df_merged_sel['ref5mer']).most_common():
        f_out.write(f'{k}\t{v}\n')

# plt.figure(figsize=(5, 5))
# plt.plot(df_merged_sel['score_miCLiP'], df_merged_sel['modRatio'], '.')
# plt.xlim([0.5, 1])
# plt.ylim([50, 100])
# plt.xlabel('miCLiP score')
# plt.ylabel('mAFiA modRatio')