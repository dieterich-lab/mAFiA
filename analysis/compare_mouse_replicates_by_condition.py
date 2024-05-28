import pandas as pd
pd.set_option('display.max_columns', 500)
from functools import reduce
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import os
from pybedtools import BedTool
from tqdm import tqdm
from collections import Counter


THRESH_CONF = 80
THRESH_COV = 20

########################################################################################################################
dataset = 'Diet'
conditions = ['WT_CD', 'WT_WD']
replicates = ['rep1', 'rep2']
# dict_condition_names = {
#     'A': 'WT_CD',
#     'B': 'M3KO_CD',
#     'C': 'WT_WD',
#     'D': 'M3KO_WD'
# }
dict_condition_colors = {
    'WT_CD': 'b',
    'M3KO_CD': 'g',
    'WT_WD': 'r',
    'M3KO_WD': 'm'
}
########################################################################################################################
# dataset = 'HFpEF'
# conditions = ['ctrl', 'HFpEF']
# replicates = ['rep1', 'rep2']
# dict_condition_colors = {
#     'ctrl': 'b',
#     'HFpEF': 'r',
# }
########################################################################################################################
# dataset = 'CM'
# conditions = ['WT', 'M3KO']
# replicates = ['rep1', 'rep2']
# dict_condition_colors = {
#     'WT': 'b',
#     'M3KO': 'r',
# }
########################################################################################################################

results_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{dataset}'
img_out = f'/home/adrian/img_out/mouse_heart/{dataset}'
os.makedirs(img_out, exist_ok=True)

gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
gtf = BedTool(gtf_file)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def import_df_results(ds_conditions, ds_replicates, thresh_confidence, thresh_coverage):
    dfs = []
    for this_cond in ds_conditions:
        for this_rep in ds_replicates:
            this_ds = f'{this_cond}_{this_rep}'
            this_df = pd.read_csv(os.path.join(results_dir, this_ds, 'chrALL.mAFiA.sites.bed'), sep='\t', dtype={'chrom': str})
            this_df = this_df[
                (this_df['confidence']>=thresh_confidence)
                * (this_df['coverage']>=thresh_coverage)
                ]
            this_df.rename(columns={
                'modRatio': f'modRatio_{this_cond}_{this_rep}',
                'coverage': f'coverage_{this_cond}_{this_rep}',
                'confidence': f'confidence_{this_cond}_{this_rep}',
            }, inplace=True)
            dfs.append(this_df)
    return reduce(lambda left, right:
                  pd.merge(left, right, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                           how='inner'), dfs)


def plot_condition_replicates(in_ax, in_df, mod_name, ds_conditions, ds_replicates, ylim=[-5, 105]):
    conds_reps = [f'{cond}_{rep}' for cond in ds_conditions for rep in ds_replicates]
    rep_labels = [rep for cond in ds_conditions for rep in ds_replicates]
    num_conds_reps = len(conds_reps)
    labelled = {cond: False for cond in ds_conditions}
    for _, this_row in in_df.iterrows():
        start = 0
        for this_cond in ds_conditions:
            this_series = [this_row[f'modRatio_{this_cond}_{this_rep}'] for this_rep in ds_replicates]
            if not labelled[this_cond]:
                ax.plot(range(start, start+len(this_series)), this_series,
                        label=this_cond,
                        c=dict_condition_colors[this_cond], linestyle='-', marker='o', alpha=0.5)
                labelled[this_cond] = True
            else:
                ax.plot(range(start, start+len(this_series)), this_series,
                        c=dict_condition_colors[this_cond], linestyle='-', marker='o', alpha=0.5)
            start += len(this_series)
        ax.plot([len(ds_replicates)-1, len(ds_replicates)],
                [this_row[f'modRatio_{ds_conditions[0]}_{ds_replicates[-1]}'], this_row[f'modRatio_{ds_conditions[1]}_{ds_replicates[0]}']],
                c='gray', linestyle='--', alpha=0.5)

        # this_series = [this_row[f'modRatio_{this_cond_rep}'] for this_cond_rep in conds_reps]
        # in_ax.plot(range(len(conds_reps)), this_series, c=color, linestyle='--', marker='o', alpha=0.5)
        # ax.plot([1, num_conditions], this_diff, c='gray', linestyle='-', alpha=0.5)
    in_ax.set_xticks(range(num_conds_reps), rep_labels)
    in_ax.set_xlim([-0.5, num_conds_reps-0.5])
    if ylim:
        ax.set_ylim(ylim)
    in_ax.set_xlabel('Sample', fontsize=12)
    in_ax.set_ylabel(f'$S_{{{dict_mod_display[mod_name]}}}$', fontsize=12)


def annotate_sites(df_in):
    base_col_order = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer', 'gene', 'region'
    ]
    df_out = []
    print('Now annotating sites...')
    for _, this_row in tqdm(df_in.iterrows()):
        annotations = gtf.intersect(BedTool.from_dataframe(this_row.to_frame().T))
        features = [this_annot.fields[2] for this_annot in annotations]
        region = [feat for feat in list(set(features)) if feat not in ['gene', 'transcript', 'exon']]
        if len(annotations):
            this_row['gene'] = annotations[0].attrs['gene_name']
            this_row['region'] = ', '.join(region)
            df_out.append(this_row)
    df_out = pd.DataFrame(df_out)
    col_order = base_col_order + \
                list(df_out.columns[df_out.columns.str.contains('modRatio')]) + \
                list(df_out.columns[df_out.columns.str.contains('coverage')]) + \
                list(df_out.columns[df_out.columns.str.contains('confidence')])
    return df_out[col_order]


def plot_bar_chart_site_trend_by_region():
    regions = ['five_prime_utr', 'CDS', 'three_prime_utr']
    # ylim = [-100, 100]
    plt.figure(figsize=(10, 5))
    all_region_site_counts = []
    for subplot_ind, mod in enumerate(['psi', 'm6A']):
        plt.subplot(1, 2, subplot_ind + 1)
        for mask_name in ['increasing', 'decreasing']:
            this_df = pd.read_csv(
                os.path.join(img_out,
                             f'{conditions[0]}_vs_{conditions[1]}_sites_{mask_name}_minCov{THRESH_COV}_deltaS{thresh_delta}.tsv'),
                sep='\t')
            all_counts = Counter(this_df[this_df['name'] == mod]['region'])
            region_site_counts = [all_counts[this_region] for this_region in regions]
            all_region_site_counts.extend(region_site_counts)
            if mask_name == 'decreasing':
                region_site_counts = [-this_count for this_count in region_site_counts]
                color = 'r'
            else:
                color = 'b'
            plt.bar(range(len(regions)), region_site_counts, width=0.5, color=color)
            plt.xticks(range(len(regions)), regions)
            plt.title(f'${dict_mod_display[mod]}$', fontsize=12)
        plt.xlabel('Transcript region', fontsize=12)
        if subplot_ind == 0:
            plt.ylabel('Sites with increaseing / decreasing trend', fontsize=12)
    ymax = np.max(all_region_site_counts)
    ytop = int(np.ceil(ymax / 25) * 25)
    for subplot_ind in range(len(['psi', 'm6A'])):
        plt.subplot(1, 2, subplot_ind + 1)
        plt.ylim([-ytop, ytop])
    plt.suptitle(f'{conditions[1]} cf. {conditions[0]}', fontsize=15)
    plt.savefig(os.path.join(img_out, f'{conditions[1]}_vs_{conditions[0]}_site_trend_by_region.png'), bbox_inches='tight')


thresh_delta = 5.0
df_merged = import_df_results(conditions, replicates, THRESH_CONF, THRESH_COV)
num_rows = len(df_merged)
modRatios = {}
for this_cond in conditions:
    modRatios[this_cond] = df_merged.loc[:, df_merged.columns.str.contains(f'modRatio_{this_cond}')].values
group_div = [this_val.max(axis=1) - this_val.min(axis=1) for this_val in modRatios.values()]
group_mean = [this_val.mean(axis=1) for this_val in modRatios.values()]
int_delta_max = np.vstack(group_div).T.max(axis=1)
# ext_delta_min = np.array([np.abs([this_row[i]-this_row[j]]).min() for this_row in np.vstack(group_mean).T for i in range(len(this_row)-1) for j in range(i+1, len(this_row))])
ext_delta = np.array([this_row[1]-this_row[0] for this_row in np.vstack(group_mean).T])

########################################################################################################################
### plot replicate differences #########################################################################################
########################################################################################################################
fig = plt.figure(figsize=(8, 8))

dict_mask = {
    'increasing': (int_delta_max < thresh_delta) * (ext_delta >= thresh_delta),
    'decreasing': (int_delta_max < thresh_delta) * (ext_delta < -thresh_delta)
}

for row_ind, mask_name in enumerate(['decreasing', 'increasing']):
    df_sel = df_merged[dict_mask[mask_name]]
    for col_ind, mod in enumerate(['m6A', 'psi']):
        subplot_ind = row_ind*2+col_ind+1
        ax = fig.add_subplot(2, 2, subplot_ind)
        plot_condition_replicates(ax, df_sel[df_sel['name']==mod], mod, conditions, replicates)
        if subplot_ind==1:
            ax.legend(loc='upper left')
    df_annot = annotate_sites(df_sel)
    df_annot.to_csv(os.path.join(img_out, f'{conditions[0]}_vs_{conditions[1]}_sites_{mask_name}_minCov{THRESH_COV}_deltaS{thresh_delta}.tsv'), sep='\t', index=False)

fig.suptitle(dataset, fontsize=15)
fig.savefig(os.path.join(img_out, f'{conditions[0]}_vs_{conditions[1]}_minCov{THRESH_COV}.png'))

plot_bar_chart_site_trend_by_region()