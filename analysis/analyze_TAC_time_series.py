import os
import pandas as pd
pd.set_option('display.max_columns', 500)
from functools import reduce
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

thresh_confidence = 50.0
thresh_coverage = 20
results_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/TAC'
img_out = '/home/adrian/img_out/TAC'
os.makedirs(img_out, exist_ok=True)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def plot_scatter(mod, group, days):
    sub_df = df_merged[df_merged['name'] == mod]
    for this_day in days:
        plt.plot(sub_df[f'modRatio_control'], sub_df[f'modRatio_{group}_day{this_day}'], '.', label=f'day{this_day}')
    plt.plot([0, 100], [0, 100], 'r--', alpha=0.5)
    plt.legend(loc='upper left')
    plt.xlabel('control', fontsize=12)
    plt.title(group, fontsize=12)
    plt.xlim([-5, 105])
    plt.ylim([-5, 105])

def plot_matchstick(ax, in_df, mod_name, group, days, color='gray', ylim=[-5, 105]):
    num_conditions = len(days)
    diffs = in_df[in_df['name'] == mod_name][[f'modRatio_{group}_{this_day}' for this_day in days]].values
    for this_diff in diffs:
        ax.plot(range(num_conditions), this_diff, c=color, marker='o', alpha=0.5, label=group)
        # ax.plot([1, num_conditions], this_diff, c='gray', linestyle='-', alpha=0.5)
    ax.set_xticks(range(num_conditions), days)
    ax.set_xlim([-0.5, num_conditions-0.5])
    if ylim:
        ax.set_ylim(ylim)
    ax.set_xlabel('condition', fontsize=12)
    ax.set_ylabel(f'$S_{{{dict_mod_display[mod_name]}}}$', fontsize=12)

dict_conditions = {
    '40-34': 'control',
    '40-29': 'TAC_day1',
    '40-26': 'TAC_day7',
    '40-33': 'TAC_day21',
    '40-30': 'TAC_day56',
    '40-31': 'SHAM_day1',
    '40-27': 'SHAM_day7',
    '40-28': 'SHAM_day21',
    '40-32': 'SHAM_day56',
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

### select sites by ts trend ###
ts_mod_ratios_SHAM = df_merged.loc[:, df_merged.columns.str.contains('modRatio_SHAM')].values
ts_mod_ratios_TAC = df_merged.loc[:, df_merged.columns.str.contains('modRatio_TAC')].values
ts_mod_ratios_diff = ts_mod_ratios_TAC - ts_mod_ratios_SHAM

# plt.figure(figsize=(10, 10))
# plt.subplot(2, 2, 1)
# plot_scatter('m6A', 'SHAM', [1, 7, 21])
# plt.subplot(2, 2, 2)
# plot_scatter('m6A', 'TAC', [1, 7, 21])
# plt.subplot(2, 2, 3)
# plot_scatter('psi', 'SHAM', [1, 7, 21])
# plt.subplot(2, 2, 4)
# plot_scatter('psi', 'TAC', [1, 7, 21])

mask_name = 'increasing'
mask = (np.diff(ts_mod_ratios_TAC, axis=1)>0).all(axis=1) * (ts_mod_ratios_diff>=20).any(axis=1)
# mask = (ts_mod_ratios_diff>0).all(axis=1) * (np.diff(ts_mod_ratios_TAC, axis=1)>0).all(axis=1)
# mask_name = 'decreasing'
# mask = (np.diff(ts_mod_ratios_TAC, axis=1)<0).all(axis=1) * (ts_mod_ratios_diff<-10).any(axis=1)
# mask_name = 'large_std'
# mask = np.std(ts_mod_ratios, axis=1)>=5.0
# mask_name = 'large_delta'
# mask = np.abs((df_merged['modRatio_day56'] - df_merged['modRatio_control']).values>=5.0)
df_sel = df_merged[mask]

# plt.figure(figsize=(10, 5))
# ax_m6A = plt.subplot(1, 2, 1)
# plot_matchstick(ax_m6A, df_sel, 'm6A', ['TAC_day1', 'TAC_day7', 'TAC_day21'])
# ax_psi = plt.subplot(1, 2, 2)
# plot_matchstick(ax_psi, df_sel, 'psi', ['TAC_day1', 'TAC_day7', 'TAC_day21'])
# plt.savefig(os.path.join(img_out, f'S_m6A_psi_time_series_{mask_name}.png'), bbox_inches='tight')
# plt.close('all')

# count = 0
for _, this_row in df_sel.iterrows():
    # count += 1
    # if count>=20:
    #     break
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot()
    plot_matchstick(ax, this_row.to_frame().T, this_row['name'], group='SHAM', days=['day1', 'day7', 'day21', 'day56'], color='blue')
    plot_matchstick(ax, this_row.to_frame().T, this_row['name'], group='TAC', days=['day1', 'day7', 'day21', 'day56'], color='red')
    ax.legend(loc='upper right')
    ax.set_title(f"chr{this_row['chrom']}: {this_row['chromEnd']}, {this_row['ref5mer'].replace('T', 'U')}")

# index = 1513
# span = 5
# sub_df = df_merged.loc[index-span:index+span]
# sub_df = sub_df[sub_df['strand']==sub_df.loc[index]['strand']]
# plt.figure(figsize=(10, 5))
# ax_m6A = plt.subplot(1, 2, 1)
# plot_matchstick(ax_m6A, sub_df, 'm6A')
# ax_psi = plt.subplot(1, 2, 2)
# plot_matchstick(ax_psi, sub_df, 'psi')

from pybedtools import BedTool
gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
gtf = BedTool(gtf_file)

col_order = [
    'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer', 'gene',
    'modRatio_control', 'modRatio_day1', 'modRatio_day7', 'modRatio_day21',	'modRatio_day56',
    'coverage_control', 'coverage_day1', 'coverage_day7', 'coverage_day21', 'coverage_day56',
    'confidence_control', 'confidence_day1', 'confidence_day7', 'confidence_day21', 'confidence_day56'
]

# annots = []
df_annot = []
for _, this_row in df_sel.iterrows():
    this_annot = gtf.intersect(BedTool.from_dataframe(this_row.to_frame().T))
    if len(this_annot):
        this_row['gene'] = this_annot[0].attrs['gene_name']
        df_annot.append(this_row)
    # annots.append(this_annot[0])
    # print(this_annot)
df_annot = pd.DataFrame(df_annot)[col_order]
df_annot.to_csv(os.path.join(img_out, f'sites_{mask_name}.tsv'), sep='\t', index=False)