import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import pysam
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
irfinder_dir = os.path.join(res_dir, 'IRFinder')

img_out = '/home/adrian/img_out/IR_junctions'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

ds = 'HFpEF'
conditions = ['ctrl_merged', 'HFpEF_merged']
# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']
# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']
cond_names = [this_cond.rstrip('_merged') for this_cond in conditions]

### delta IR ratio ###
# this_cond = conditions[0]
cond_df_irf = {}
for this_cond in conditions:
    df_irf = pd.read_csv(os.path.join(irfinder_dir, f'{ds}_{this_cond}.txt'), sep='\t', dtype={'Chr': str})
    df_irf = df_irf.rename(columns={
        'Chr': 'chrom',
        'Start': 'chromStart',
        'End': 'chromEnd',
        'Name': 'name',
        'Null': 'score',
        'Strand': 'strand'
    })
    df_irf_clean = df_irf.loc[['known-exon' not in this_name.split('/')[-1] for this_name in df_irf['name']], :]
    # cond_df_irf[this_cond] = df_irf_clean[df_irf_clean['Warnings'] == '-']
    cond_df_irf[this_cond] = df_irf_clean

thresh_IRratio = 0.01
thresh_IntronDepth = 5

df_irf_merged = pd.merge(*list(cond_df_irf.values()), on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'],
                         suffixes=['_'+this_cond_name for this_cond_name in cond_names])
# df_irf_merged_thresh = df_irf_merged[
#     ((df_irf_merged.loc[:, df_irf_merged.keys().str.find('IRratio') >= 0] >= thresh_IRratio).any(axis=1))
#     * ((df_irf_merged.loc[:, df_irf_merged.keys().str.find('IntronDepth_') >= 0] >= thresh_IntronDepth).any(axis=1))
# ].copy()
df_irf_merged_thresh = df_irf_merged[
    # ((df_irf_merged.loc[:, df_irf_merged.keys().str.find('IntronDepth_') >= 0] >= thresh_IntronDepth).any(axis=1))
    (df_irf_merged[f'IntronDepth_{cond_names[0]}'] >= thresh_IntronDepth)
    * (df_irf_merged[f'IntronDepth_{cond_names[1]}'] >= thresh_IntronDepth)
    * (df_irf_merged[f'IRratio_{cond_names[0]}'] >= thresh_IRratio)
    * (df_irf_merged[f'IRratio_{cond_names[0]}'] < 1.0)
    * (df_irf_merged[f'IRratio_{cond_names[1]}'] >= thresh_IRratio)
    * (df_irf_merged[f'IRratio_{cond_names[1]}'] < 1.0)
    ].copy()
# df_irf_merged_thresh['delta_IRratio'] = df_irf_merged_thresh[f'IRratio_{cond_names[1]}'] - \
#                                         df_irf_merged_thresh[f'IRratio_{cond_names[0]}']
# df_irf_merged_thresh.sort_values('delta_IRratio', inplace=True)
df_irf_merged_thresh['log2fc_IRratio'] = np.log2(df_irf_merged_thresh[f'IRratio_{cond_names[1]}'] /
                                                 df_irf_merged_thresh[f'IRratio_{cond_names[0]}'])
df_irf_merged_thresh.sort_values('log2fc_IRratio', inplace=True)


# df_irf_merged_thresh.to_csv(
#     os.path.join(img_out, f"ir_list_{'_'.join(cond_names)}_IRratio{thresh_IRratio}_IntronDepth{thresh_IntronDepth}.tsv"),
#     sep='\t',
#     index=False
# )

### scatter plot of IR ratios ###
xylim = [-0.05, 1.05]
plt.figure(figsize=(5, 5))
plt.plot(df_irf_merged_thresh[f'IRratio_{cond_names[0]}'], df_irf_merged_thresh[f'IRratio_{cond_names[1]}'], '.')
plt.plot([0, 1], [0, 1], c='gray', ls='--')
plt.xlim(xylim)
plt.ylim(xylim)
plt.xlabel(f'IR ratio, {cond_names[0]}')
plt.ylabel(f'IR ratio, {cond_names[1]}')

### delta mod ###
# thresh_conf = 0.0
# cond_df_mod = {}
# for this_cond in conditions:
#     cond_df_mod[this_cond] = pd.read_csv(os.path.join(res_dir, ds, this_cond, f'chrALL.mAFiA.sites.bed'),
#                                          sep='\t', dtype={'chrom': str})
# df_mod_merged = pd.merge(*list(cond_df_mod.values()),  on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'],
#                          suffixes=['_'+this_cond_name for this_cond_name in cond_names])
# df_mod_merged_thresh = df_mod_merged[
#     (df_mod_merged[f'confidence_{cond_names[0]}'] >= thresh_conf)
#     * (df_mod_merged[f'confidence_{cond_names[1]}'] >= thresh_conf)
# ].copy()
# df_mod_merged_thresh['delta_modRatio'] = df_mod_merged_thresh[f'modRatio_{cond_names[1]}'] - df_mod_merged_thresh[f'modRatio_{cond_names[0]}']

def get_margin_avg_mod_prob(in_bam, in_row, margin=50, thresh_coverage=5):
    flag_required = 0 if in_row['strand'] == '+' else 16
    mod_read_mod_probs = {this_mod: [] for this_mod in mods}
    read_count = 0
    for this_read in in_bam.fetch(in_row['chrom'],
                                  in_row['chromStart'] - margin,
                                  in_row['chromEnd'] + margin):
        if this_read.flag == flag_required:
            query_to_ref_pos = {pair[0]: pair[1] for pair in this_read.get_aligned_pairs(matches_only=True)}
            for mod in mods:
                ref_pos_mod_prob = [(query_to_ref_pos[pair[0]], pair[1])
                                    for pair in this_read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])
                                    ]
                left_ref_pos_mod_prob = [pair for pair in ref_pos_mod_prob if
                                         (pair[0] >= (in_row['chromStart'] - margin))
                                         and (pair[0] <= in_row['chromStart'])]
                right_ref_pos_mod_prob = [pair for pair in ref_pos_mod_prob if
                                          (pair[0] >= in_row['chromEnd'])
                                          and (pair[0] <= (in_row['chromEnd'] + margin))]
                if len(left_ref_pos_mod_prob) and len(right_ref_pos_mod_prob):
                    read_count += 1
                    mod_read_mod_probs[mod].extend(left_ref_pos_mod_prob + right_ref_pos_mod_prob)
    # print(read_count)
    if read_count >= thresh_coverage:
        return {mod: np.mean([pair[1] >= 128 for pair in mod_read_mod_probs[mod]]) for mod in mods}
    else:
        return {mod: np.nan for mod in mods}

index_cond_avg_mod_prob = {this_index: {} for this_index in df_irf_merged_thresh.index}
for this_cond in conditions:
    this_cond_bam_path = os.path.join(res_dir, ds, this_cond, 'chrALL.mAFiA.reads.bam')
    with pysam.AlignmentFile(this_cond_bam_path, 'rb') as this_cond_bam:
        for this_index, this_row in tqdm(df_irf_merged_thresh.iterrows()):
            index_cond_avg_mod_prob[this_index][this_cond] = get_margin_avg_mod_prob(this_cond_bam, this_row)

for this_mod in mods:
    # df_irf_merged_thresh[f'delta_{this_mod}'] = [v[conditions[1]][this_mod] - v[conditions[0]][this_mod]
    #                                              for k, v in index_cond_avg_mod_prob.items()]
    df_irf_merged_thresh[f'log2fc_{this_mod}'] = [np.log2(v[conditions[1]][this_mod] / v[conditions[0]][this_mod])
                                                 for k, v in index_cond_avg_mod_prob.items()]

print(df_irf_merged_thresh[['log2fc_IRratio', 'log2fc_m6A', 'log2fc_psi']])

# xmax = np.round((np.nanmax(np.abs(df_irf_merged_thresh['delta_IRratio'].values)) // 0.05 + 1) * 0.05, 2)
# ymax = np.round((np.nanmax(np.abs(df_irf_merged_thresh[['log2fc_m6A', 'log2fc_psi']].values)) // 0.05 + 1) * 0.05, 2)
xmax = 1.0
ymax = 0.8

# xticks = np.round(np.arange(-xmax, xmax+0.01, 0.05), 2)
# yticks = np.round(np.arange(-ymax, ymax+0.01, 0.05), 2)
xticks = np.round(np.linspace(-xmax, xmax, 5), 2)
yticks = np.round(np.linspace(-ymax, ymax, 5), 2)

bin_centers = 0.5 * (xticks[1:] + xticks[:-1])

boundary = 0.1

global_IR_shift = np.median([x for x in df_irf_merged_thresh[f'log2fc_IRratio'].values if ~np.isnan(x) and ~np.isinf(x)])

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    # plt.plot(df_irf_merged_thresh[f'delta_{this_mod}'], df_irf_merged_thresh['delta_IRratio'], 'o')
    # plt.axvline(x=0, c='gray')
    # plt.axhline(y=0, c='gray')
    binned_y = []
    vec_y = df_irf_merged_thresh[f'log2fc_{this_mod}'].values
    vec_x = df_irf_merged_thresh[f'log2fc_IRratio'].values
    # vec_x = df_irf_merged_thresh['log2fc_IRratio'].values
    # for this_bin_ind in range(len(xticks)-1):
    #     bin_start = xticks[this_bin_ind]
    #     bin_end = xticks[this_bin_ind+1]
    #     mask = (vec_x >= bin_start) * (vec_x < bin_end)
    #     binned_y.append([x for x in df_irf_merged_thresh[f'log2fc_{this_mod}'].values[mask] if ~np.isnan(x)])
    # binned_y = [
    #     [x for x in df_irf_merged_thresh[f'log2fc_{this_mod}'].values[vec_x < -boundary] if ~np.isnan(x) and ~np.isinf(x)],
    #     [x for x in df_irf_merged_thresh[f'log2fc_{this_mod}'].values[vec_x > boundary] if ~np.isnan(x) and ~np.isinf(x)],
    # ]
    binned_y = [
        # [x for x in df_irf_merged_thresh[f'log2fc_IRratio'].values[(vec_x < -boundary) * (vec_x >= -1.0)] if ~np.isnan(x) and ~np.isinf(x)],
        # [x for x in df_irf_merged_thresh[f'log2fc_IRratio'].values[(vec_x >= boundary) * (vec_x < 1.0)] if ~np.isnan(x) and ~np.isinf(x)],
        [x for x in vec_y[(vec_x < -boundary)] if
         ~np.isnan(x) and ~np.isinf(x)],
        [x for x in vec_y[(vec_x >= boundary)] if
         ~np.isnan(x) and ~np.isinf(x)],
    ]
    print(this_mod, [len(ll) for ll in binned_y])
    # plt.axhline(y=global_IR_shift, c='red', ls='--', alpha=0.5)
    plt.boxplot(binned_y, positions=[-0.5, 0.5], widths=0.25, whis=1.5, showfliers=False)
    # plt.violinplot(binned_y, positions=[-0.5, 0.5])
    # plt.bar([-0.5, 0.5], [np.median(this_bin) for this_bin in binned_y], width=0.3)
    plt.xlim([-xmax, xmax])
    # plt.ylim([-ymax, ymax])
    # plt.xticks(xticks, xticks)
    plt.xticks([-0.5, 0.5], [f'< -{boundary} ({len(binned_y[0])})', rf'$\geq$ {boundary} ({len(binned_y[1])})'])
    plt.yticks(yticks, yticks)
    plt.ylabel(r'$log_{2}fc$ $\bar{S}$')
    plt.xlabel(r'$log_{2}fc$ IR ratio')
    plt.title(rf'${{{dict_mod_display[this_mod]}}}$')
plt.suptitle(f'{ds}\n{cond_names[1]} vs {cond_names[0]}')
# plt.savefig(os.path.join(img_out, f'log2fc_IR_vs_log2fc_mod_{ds}.png'), bbox_inches='tight')
plt.savefig(os.path.join(img_out, f'log2fc_mod_vs_log2fc_IR_{ds}.png'), bbox_inches='tight')