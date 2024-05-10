import os
import pandas as pd
pd.set_option('display.max_columns', 500)
from functools import reduce
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

thresh_confidence = 80.0
thresh_coverage = 20
results_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/TAC'
img_out = '/home/adrian/img_out/TAC'
os.makedirs(img_out, exist_ok=True)

dict_conditions = {
    '40-34': 'control',
    '40-29': 'day1',
    '40-26': 'day7',
    # '40-33': 'day21',
    # '40-30': 'day56'
}
conditions = list(dict_conditions.values())

dfs = []
for this_ds in dict_conditions.keys():
    this_df = pd.read_csv(os.path.join(results_dir, this_ds, 'chrALL.mAFiA.sites.bed'), sep='\t', dtype={'chrom': str})
    this_df = this_df[
        (this_df['confidence']>=thresh_confidence)
        * (this_df['coverage']>=thresh_coverage)
        ]
    this_df.rename(columns={
        'modRatio': f'modRatio_{dict_conditions[this_ds]}',
        'coverage': f'coverage_{dict_conditions[this_ds]}',
        'confidence': f'confidence_{dict_conditions[this_ds]}',
    }, inplace=True)
    dfs.append(this_df)
df_merged = reduce(lambda left, right:
                   pd.merge(left, right, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                            how='inner'), dfs)
df_merged['delta'] = df_merged['modRatio_day7'] - df_merged['modRatio_control']
df_large_delta = df_merged[np.abs(df_merged['delta'].values)>=10.0]

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def plot_matchstick(ax, in_df, mod_name, ylim=[-5, 105]):
    num_conditions = len(conditions)
    diffs = in_df[in_df['name'] == mod_name][[f'modRatio_{this_cond}' for this_cond in conditions]].values
    for this_diff in diffs:
        ax.plot(range(num_conditions), this_diff, c='gray', marker='o', alpha=0.5)
        # ax.plot([1, num_conditions], this_diff, c='gray', linestyle='-', alpha=0.5)
    ax.set_xticks(range(num_conditions), conditions)
    ax.set_xlim([-0.5, num_conditions-0.5])
    ax.set_ylim(ylim)
    ax.set_xlabel('condition', fontsize=12)
    ax.set_ylabel(f'$S_{{{dict_mod_display[mod_name]}}}$', fontsize=12)

plt.figure(figsize=(10, 5))
ax_m6A = plt.subplot(1, 2, 1)
plot_matchstick(ax_m6A, df_large_delta, 'm6A')
ax_psi = plt.subplot(1, 2, 2)
plot_matchstick(ax_psi, df_large_delta, 'psi')
plt.savefig(os.path.join(img_out, 'S_m6A_psi_time_series.png'), bbox_inches='tight')
plt.close('all')