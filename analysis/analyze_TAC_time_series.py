import os
import pandas as pd
pd.set_option('display.max_columns', 500)
from functools import reduce
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pybedtools import BedTool
from collections import Counter


thresh_confidence = 50.0
thresh_coverage = 20
results_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v0/TAC'
img_out = '/home/adrian/img_out/TAC'
os.makedirs(img_out, exist_ok=True)

gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
gtf = BedTool(gtf_file)

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
    # '40-34': 'control',
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
mask = (np.diff(ts_mod_ratios_TAC, axis=1)>0).all(axis=1)\
       * (ts_mod_ratios_diff>=10).any(axis=1)\
       * ((ts_mod_ratios_SHAM.max(axis=1) - ts_mod_ratios_SHAM.min(axis=1)) < 10)

# mask_name = 'decreasing'
# mask = (np.diff(ts_mod_ratios_TAC, axis=1)<0).all(axis=1)\
#        * (ts_mod_ratios_diff<-10).any(axis=1) \
#        * ((ts_mod_ratios_SHAM.max(axis=1) - ts_mod_ratios_SHAM.min(axis=1)) < 10)

df_sel = df_merged[mask]

num_display = 16
num_rows = 4
num_cols = 4
if mask_name == 'increasing':
    display_indices = (df_sel['modRatio_TAC_day56'] - df_sel['modRatio_SHAM_day56']).nlargest(num_display).index
elif mask_name == 'decreasing':
    display_indices = (df_sel['modRatio_SHAM_day56'] - df_sel['modRatio_TAC_day56']).nlargest(num_display).index

fig = plt.figure(figsize=(15, 15))
for subplot_ind, this_ind in enumerate(display_indices):
    this_row = df_sel.loc[this_ind]
    ax = fig.add_subplot(num_rows, num_cols, subplot_ind+1)
    plot_matchstick(ax, this_row.to_frame().T, this_row['name'], group='SHAM', days=['day1', 'day7', 'day21', 'day56'], color='blue')
    plot_matchstick(ax, this_row.to_frame().T, this_row['name'], group='TAC', days=['day1', 'day7', 'day21', 'day56'], color='red')
    ax.legend(loc='upper left', title=f"chr{this_row['chrom']}: {this_row['chromEnd']}\n{this_row['ref5mer'].replace('T', 'U')}")
    # ax.set_title(f"chr{this_row['chrom']}: {this_row['chromEnd']}, {this_row['ref5mer'].replace('T', 'U')}")
fig.suptitle(f'TAC {mask_name}', fontsize=15)
fig.savefig(os.path.join(img_out, f'ts_{mask_name}.png'), bbox_inches='tight')

### annotate sites in gene region ###
col_order = [
    'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer', 'gene', 'region',
    'modRatio_TAC_day1', 'modRatio_TAC_day7', 'modRatio_TAC_day21',	'modRatio_TAC_day56',
    'modRatio_SHAM_day1', 'modRatio_SHAM_day7', 'modRatio_SHAM_day21', 'modRatio_SHAM_day56',
    'coverage_TAC_day1', 'coverage_TAC_day7', 'coverage_TAC_day21', 'coverage_TAC_day56',
    'coverage_SHAM_day1', 'coverage_SHAM_day7', 'coverage_SHAM_day21', 'coverage_SHAM_day56',
    'confidence_TAC_day1', 'confidence_TAC_day7', 'confidence_TAC_day21', 'confidence_TAC_day56',
    'confidence_SHAM_day1', 'confidence_SHAM_day7', 'confidence_SHAM_day21', 'confidence_SHAM_day56'
]

df_annot = []
for _, this_row in df_sel.iterrows():
    annotations = gtf.intersect(BedTool.from_dataframe(this_row.to_frame().T))
    features = [this_annot.fields[2] for this_annot in annotations]
    region = [feat for feat in list(set(features)) if feat not in ['gene', 'transcript', 'exon']]
    # print(region)
    if len(annotations):
        this_row['gene'] = annotations[0].attrs['gene_name']
        this_row['region'] = ', '.join(region)
        df_annot.append(this_row)
df_annot = pd.DataFrame(df_annot)[col_order]
df_annot.to_csv(os.path.join(img_out, f'sites_{mask_name}.tsv'), sep='\t', index=False)

### plot sites with increasing / decreasing trend ###
regions = ['five_prime_utr', 'CDS', 'three_prime_utr']
ylim = [-60, 60]
plt.figure(figsize=(10, 4))
for subplot_ind, mod in enumerate(['psi', 'm6A']):
    plt.subplot(1, 2, subplot_ind+1)
    for mask_name in ['increasing', 'decreasing']:
        this_df = pd.read_csv(os.path.join(img_out, f'sites_{mask_name}.tsv'), sep='\t')
        all_counts = Counter(this_df[this_df['name']==mod]['region'])
        region_site_counts = [all_counts[this_region] for this_region in regions]
        if mask_name=='decreasing':
            region_site_counts = [-this_count for this_count in region_site_counts]
            color = 'r'
        else:
            color = 'b'
        plt.bar(range(len(regions)), region_site_counts, width=0.5, color=color)
        plt.xticks(range(len(regions)), regions)
        plt.title(f'${dict_mod_display[mod]}$', fontsize=12)
    plt.xlabel('Transcript region', fontsize=12)
    if subplot_ind==0:
        plt.ylabel('Sites with increaseing / decreasing trend', fontsize=12)
    plt.ylim(ylim)
plt.savefig(os.path.join(img_out, f'barchart_sites_with_increasing_decreasing_trend.png'), bbox_inches='tight')