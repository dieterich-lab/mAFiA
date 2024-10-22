import os
import pandas as pd
import numpy as np
import matplotlib as mpl
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
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################
import matplotlib.pyplot as plt

bed6_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

miCLIP_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/miCLIP'
mAFiA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'
os.makedirs(img_out, exist_ok=True)

# cond = 'SHAM'
# cond = 'TAC'

conditions = ['SHAM', 'TAC']

for cond in conditions:

    miCLIP_file = os.path.join(miCLIP_dir, f'nochr_predicted_m6A_{cond}.bed')
    mAFiA_file = os.path.join(mAFiA_dir, f'{cond}_merged', 'chrALL.mAFiA.sites.bed')

    df_miCLIP = pd.read_csv(miCLIP_file, sep='\t', names=bed6_fields)
    df_miCLIP['name'] = 'm6A'
    df_mAFiA = pd.read_csv(mAFiA_file, sep='\t', dtype={'chrom': str})
    df_mAFiA.drop(columns=['score'], inplace=True)

    df_merged = pd.merge(df_miCLIP, df_mAFiA, on=[
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'strand'
    ])

    print(f'Enough coverage for {len(df_merged)} / {len(df_miCLIP)} miCLIP sites')

    thresh_confidence = 80.0
    df_merged_thresh = df_merged[df_merged['confidence'] >= thresh_confidence]

    # plt.figure(figsize=(5, 5))
    # plt.scatter(df_merged_thresh['score'], df_merged_thresh['modRatio'], s=1)

    bin_min = 0.5
    bin_max = 1.0
    bin_width = 0.1
    bin_edges = np.arange(bin_min, bin_max+0.01, bin_width)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    binned_modRatios = []
    for bin_ind in range(len(bin_edges)-1):
        bin_start = bin_edges[bin_ind]
        bin_end = bin_edges[bin_ind+1] + 0.000001
        binned_modRatios.append(
            df_merged_thresh[(df_merged_thresh['score'] >= bin_start) * (df_merged_thresh['score'] < bin_end)]['modRatio'].values)

    plt.figure(figsize=(5, 5))
    plt.boxplot(binned_modRatios, positions=bin_centers, widths=0.05, whis=0.5, showfliers=False)
    plt.xlim([bin_min, bin_max])
    plt.xticks(bin_edges, np.round(bin_edges, 1))
    # plt.xlabel('miCLIP score', fontsize=12)
    # plt.ylabel('mAFiA stoichiometry', fontsize=12)
    # plt.title(cond, fontsize=15)
    plt.savefig(os.path.join(img_out, f'boxplot_mAFiA_miCLIP_{cond}.{FMT}'), **fig_kwargs)

### compare miCLIP diff. sites to mAFiA KS sites ###
from matplotlib_venn import venn2

# thresh_pval = 0.05
# ks_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/KS_test/TAC_SHAM_merged_TAC_merged.bed'
# df_ks = pd.read_csv(ks_file, sep='\t', dtype={'chrom': str})
# df_ks_thresh = df_ks[df_ks['pval'] < thresh_pval]

# df_mAFiA_sham = pd.read_csv(os.path.join(mAFiA_dir, f'SHAM_merged', 'chrALL.mAFiA.sites.bed'),
#                             sep='\t', dtype={'chrom': str})
# df_mAFiA_tac = pd.read_csv(os.path.join(mAFiA_dir, f'TAC_merged', 'chrALL.mAFiA.sites.bed'),
#                             sep='\t', dtype={'chrom': str})
#
# df_miCLIP_sham = pd.read_csv(os.path.join(miCLIP_dir, f'nochr_predicted_m6A_SHAM.bed'), sep='\t', names=bed6_fields)
# df_miCLIP_sham['name'] = 'm6A'
# df_miCLIP_tac = pd.read_csv(os.path.join(miCLIP_dir, f'nochr_predicted_m6A_TAC.bed'), sep='\t', names=bed6_fields)
# df_miCLIP_tac['name'] = 'm6A'
#
# df_miCLIP_merged = pd.merge(df_miCLIP_sham, df_miCLIP_tac, on=[
#     'chrom',
#     'chromStart',
#     'chromEnd',
#     'name',
#     'strand'
# ], suffixes=['_SHAM', '_TAC'])
#
# thresh = 0.9
# miCLIP_plus_sites = [tuple(val) for val in df_miCLIP_merged[
#     (df_miCLIP_merged['score_SHAM'] < thresh)
#     * (df_miCLIP_merged['score_TAC'] >= thresh)
#     * (df_miCLIP_merged['score_delta'] >= 0.1)
#     ][['chrom', 'chromStart']].values]
# modRatios_sham_tac = []
# for this_site in miCLIP_plus_sites:
#     sub_df_mAFiA_sham = df_mAFiA_sham[
#         (df_mAFiA_sham['chrom'] == this_site[0])
#         * (df_mAFiA_sham['chromStart'] == this_site[1])
#     ]
#     sub_df_mAFiA_tac = df_mAFiA_tac[
#         (df_mAFiA_tac['chrom'] == this_site[0])
#         * (df_mAFiA_tac['chromStart'] == this_site[1])
#     ]
#     if len(sub_df_mAFiA_sham) and len(sub_df_mAFiA_tac):
#         modRatio_sham = sub_df_mAFiA_sham['modRatio'].values[0]
#         modRatio_tac = sub_df_mAFiA_tac['modRatio'].values[0]
#         modRatios_sham_tac.append((modRatio_sham, modRatio_tac))


# df_miCLIP_merged['score_delta'] = df_miCLIP_merged['score_TAC'] - df_miCLIP_merged['score_SHAM']
# df_miCLIP_merged_ks = pd.merge(df_miCLIP_merged, df_ks, on=[
#     'chrom',
#     'chromStart',
#     'chromEnd',
#     'name',
#     'strand'
# ])
#
# miCLIP_diff_sites_pos = set([
#     tuple(this_val)
#     for this_val in df_miCLIP_merged_ks[df_miCLIP_merged_ks['score_delta'] > 0.0][['chrom', 'chromStart']].values
# ])
#
# mAFiA_diff_sites_pos = set([
#     tuple(this_val)
#     for this_val in df_miCLIP_merged_ks[
#         (df_miCLIP_merged_ks['delta'] > 0.0)
#         ][['chrom', 'chromStart']].values
# ])
#
# miCLIP_diff_sites_neg = set([
#     tuple(this_val)
#     for this_val in df_miCLIP_merged_ks[df_miCLIP_merged_ks['score_delta'] < 0.0][['chrom', 'chromStart']].values
# ])
#
# mAFiA_diff_sites_neg = set([
#     tuple(this_val)
#     for this_val in df_miCLIP_merged_ks[
#         (df_miCLIP_merged_ks['delta'] < 0.0)
#         ][['chrom', 'chromStart']].values
# ])
#
# plt.figure(figsize=(5, 5))
# plt.subplot(1, 2, 1)
# venn2([miCLIP_diff_sites_pos, mAFiA_diff_sites_pos], ['miCLIP', 'mAFiA'], set_colors=['blue', 'red'])
# plt.title('+ve $\Delta$', fontsize=15)
# plt.subplot(1, 2, 2)
# venn2([miCLIP_diff_sites_neg, mAFiA_diff_sites_neg], ['miCLIP', 'mAFiA'], set_colors=['blue', 'red'])
# plt.title('-ve $\Delta$', fontsize=15)

# plt.figure(figsize=(5, 5))
# plt.scatter(df_miCLIP_merged_ks['score_delta'], df_miCLIP_merged_ks['delta'], s=2)
#
# bin_min = -0.3
# bin_max = 0.3
# bin_width = 0.1
# bin_edges = np.round(np.arange(bin_min, bin_max+0.01, bin_width), 1)
# bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
#
# binned_deltas = []
# for bin_ind in range(len(bin_edges)-1):
#     bin_start = bin_edges[bin_ind]
#     bin_end = bin_edges[bin_ind+1] + 0.000001
#     binned_deltas.append(
#         df_miCLIP_merged_ks[
#             (df_miCLIP_merged_ks['score_delta'] >= bin_start)
#             * (df_miCLIP_merged_ks['score_delta'] < bin_end)]
#         ['delta'].values)
#
# plt.figure(figsize=(5, 5))
# plt.boxplot(binned_deltas, positions=bin_centers, widths=0.05, showfliers=False)
# plt.xlim([bin_min, bin_max])
# plt.xticks(bin_edges, np.round(bin_edges, 1))
# plt.xlabel('change miCLIP score', fontsize=12)
# plt.ylabel('change mAFiA stoichiometry', fontsize=12)
