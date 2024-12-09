import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import pysam
import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200, transparent=True)
#######################################################################


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

def get_margin_avg_logit(in_bam, in_row, margin=50, thresh_coverage=1):
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
        return {mod: np.mean([np.log2(pair[1] / (256 - pair[1])) for pair in mod_read_mod_probs[mod]]) for mod in mods}
    else:
        return {mod: np.nan for mod in mods}


def get_logit(in_val):
    return np.log2(in_val / (1.0 - in_val))

def get_df_irf_merged_thresh():
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

    df_irf_merged = pd.merge(*list(cond_df_irf.values()),
                             on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'],
                             suffixes=['_' + this_cond_name for this_cond_name in cond_names])
    df_out = df_irf_merged[
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
    df_out['log2fc_IRratio'] = np.log2(df_out[f'IRratio_{cond_names[1]}'] /
                                                     df_out[f'IRratio_{cond_names[0]}'])
    df_out.sort_values('log2fc_IRratio', inplace=True)

    index_cond_avg_mod_prob = {this_index: {} for this_index in df_out.index}
    for this_cond in conditions:
        this_cond_bam_path = os.path.join(res_dir, ds, this_cond, 'chrALL.mAFiA.reads.bam')
        with pysam.AlignmentFile(this_cond_bam_path, 'rb') as this_cond_bam:
            for this_index, this_row in tqdm(df_out.iterrows()):
                index_cond_avg_mod_prob[this_index][this_cond] = get_margin_avg_mod_prob(this_cond_bam, this_row)
                # index_cond_avg_mod_prob[this_index][this_cond] = get_margin_avg_logit(this_cond_bam, this_row)

    # for mod in mods:
    #     df_out[f'log2fc_{mod}'] = [np.log2(v[conditions[1]][mod] / v[conditions[0]][mod])
    #                                                   for k, v in index_cond_avg_mod_prob.items()]
    #
    # print(df_out[['log2fc_IRratio', 'log2fc_m6A', 'log2fc_psi']])
    # return df_out

    for mod in mods:
        df_out[f'delta_logit_{mod}'] = [get_logit(v[conditions[1]][mod]) - get_logit(v[conditions[0]][mod])
                                        for k, v in index_cond_avg_mod_prob.items()]

    print(df_out[['log2fc_IRratio', 'delta_logit_m6A', 'delta_logit_psi']])
    return df_out

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
irfinder_dir = os.path.join(res_dir, 'IRFinder')

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'
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

# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']
# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']
ds = 'Diet'
conditions = ['WT_CD_merged', 'WT_WD_merged']

cond_names = [this_cond.rstrip('_merged') for this_cond in conditions]

thresh_IRratio = 0.01
thresh_IntronDepth = 5

df_irf_merged_thresh_filename = os.path.join(
    irfinder_dir,
    f'df_irf_merged_{ds}_IRratio{thresh_IRratio}_IntronDepth{thresh_IntronDepth}.tsv'
)
# if os.path.exists(df_irf_merged_thresh_filename):
#     df_irf_merged_thresh = pd.read_csv(df_irf_merged_thresh_filename, sep='\t')
# else:
#     df_irf_merged_thresh = get_df_irf_merged_thresh()
#     df_irf_merged_thresh.to_csv(df_irf_merged_thresh_filename, sep='\t', index=False, float_format='%.6f')
df_irf_merged_thresh = get_df_irf_merged_thresh()

### scatter plot of IR ratios ###
df_irf_merged_thresh['delta_IRratio'] = df_irf_merged_thresh[f'IRratio_{cond_names[1]}'] - df_irf_merged_thresh[f'IRratio_{cond_names[0]}']
thresh_delta_IRratio = 0.1
mask_positive = df_irf_merged_thresh['delta_IRratio'].values >= thresh_delta_IRratio
mask_negative = df_irf_merged_thresh['delta_IRratio'].values < -thresh_delta_IRratio
up_reg_genes = [this_name.split('/')[0] for this_name in df_irf_merged_thresh.loc[mask_positive, 'name'].values]
down_reg_genes = [this_name.split('/')[0] for this_name in df_irf_merged_thresh.loc[mask_negative, 'name'].values]

xylim = [0, 1.0]
xyticks = np.linspace(0, 1, 5)

plt.figure(figsize=(4*cm, 4*cm))
plt.scatter(df_irf_merged_thresh[f'IRratio_{cond_names[0]}'], df_irf_merged_thresh[f'IRratio_{cond_names[1]}'], s=1, c='gray', alpha=0.5)
plt.scatter(df_irf_merged_thresh.loc[mask_positive, f'IRratio_{cond_names[0]}'],
            df_irf_merged_thresh.loc[mask_positive, f'IRratio_{cond_names[1]}'],
            s=1, fc='r')
plt.scatter(df_irf_merged_thresh.loc[mask_negative, f'IRratio_{cond_names[0]}'],
            df_irf_merged_thresh.loc[mask_negative, f'IRratio_{cond_names[1]}'],
            s=1, fc='b')
plt.plot([0, 1], [0, 1], c='gray', ls='--')
plt.text(0.01, 0.99, '\n'.join(up_reg_genes), va='top', ha='left', c='red')
plt.text(0.99, 0.01, '\n'.join(down_reg_genes), va='bottom', ha='right', c='blue')
plt.xlim(xylim)
plt.ylim(xylim)
plt.xticks(xyticks)
plt.yticks(xyticks)
# plt.xlabel(f'IR ratio, {cond_names[0]}')
# plt.ylabel(f'IR ratio, {cond_names[1]}')
plt.savefig(os.path.join(img_out, f'scatterplot_log2fc_IRratio_{cond_names[1]}_vs_{cond_names[0]}.{FMT}'), **fig_kwargs)

### boxplot of log2fc stoichiometry versus log2fc IR ratios ###
xmax = 1.0
ymax = 0.8
xticks = np.round(np.linspace(-xmax, xmax, 5), 2)
yticks = np.round(np.linspace(-ymax, ymax, 5), 2)
bin_centers = 0.5 * (xticks[1:] + xticks[:-1])
boundary = 0.1

plt.figure(figsize=(5*cm, 4*cm))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    binned_y = []
    # vec_y = df_irf_merged_thresh[f'log2fc_{this_mod}'].values
    vec_y = df_irf_merged_thresh[f'delta_logit_{this_mod}'].values
    vec_x = df_irf_merged_thresh[f'log2fc_IRratio'].values
    binned_y = [
        [x for x in vec_y[(vec_x < -boundary)] if
         ~np.isnan(x) and ~np.isinf(x)],
        [x for x in vec_y[(vec_x >= boundary)] if
         ~np.isnan(x) and ~np.isinf(x)],
    ]
    plt.boxplot(binned_y, positions=[-0.5, 0.5], widths=0.25, showfliers=False)
    plt.xlim([-xmax, xmax])
    # plt.xticks([-0.5, 0.5], [f'< -{boundary} ({len(binned_y[0])})', rf'$\geq$ {boundary} ({len(binned_y[1])})'])
    plt.xticks([-0.5, 0.5], [f'< -{boundary}', rf'$\geq$ {boundary}'])
    if mod_ind == 0:
        plt.yticks(yticks, yticks)
    else:
        plt.yticks(yticks, [])
    # plt.ylabel(r'$log_{2}fc$ $\bar{S}$')
    # plt.xlabel(r'$log_{2}fc$ IR ratio')
    # plt.title(rf'${{{dict_mod_display[this_mod]}}}$')
# plt.suptitle(f'{ds}\n{cond_names[1]} vs {cond_names[0]}')
plt.savefig(os.path.join(img_out, f'boxplot_delta_logit_stoichiometry_vs_log2fc_IRratio_{ds}.{FMT}'), **fig_kwargs)