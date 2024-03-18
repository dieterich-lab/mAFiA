import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
# FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/HEK293/100_WT_0_IVT'
# chrs = [str(i) for i in range(2, 11)] + ['13', '14', '18', '20']
# chrs = list(range(1, 23)) + ['X']
chrs = [str(i) for i in range(2, 8)]
img_out = '/home/adrian/NCOMMS_revision/images/GLORI_all_pval'
os.makedirs(img_out, exist_ok=True)

glori_path = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.ref5mers.csv'
miclip_path = '/home/adrian/Data/DRACH/miCLIP_union_flat_exclude_Y_chromosome.ref5mers.bed'

# blacklist_motifs = ['AGACC', 'AAACA', 'TAACA', 'TAACT', 'GAACA', 'TGACC']
blacklist_motifs = []

def import_mAFiA(thresh_coverage=10):
    # dfs = []
    # for this_chr in chrs:
    #     dfs.append(pd.read_csv(os.path.join(source_data_dir, f'chr{this_chr}', 'mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t'))
    # df_mAFiA = pd.concat(dfs, ignore_index=True)
    df_mAFiA = pd.read_csv(os.path.join(source_data_dir, 'merged.mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t')
    df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage']>=thresh_coverage]
    df_mAFiA_thresh = df_mAFiA_thresh[~df_mAFiA_thresh['ref5mer'].isin(blacklist_motifs)]
    return df_mAFiA_thresh


def import_glori(thresh_pval=1.0):
    df_glori = pd.read_csv(glori_path, dtype={'Chr': str})
    df_glori['chrom'] = [chr.lstrip('chr') for chr in df_glori['Chr']]
    df_glori['chromStart'] = df_glori['Sites'] - 1
    df_glori['modRatio'] = np.int32(np.round(df_glori['NormeRatio'] * 100.0))
    df_glori_thresh = df_glori[df_glori['P_adjust']<thresh_pval]
    df_glori_thresh = df_glori_thresh[~df_glori_thresh['ref5mer'].isin(blacklist_motifs)]

    return df_glori_thresh


df_mAFiA = import_mAFiA()
df_glori = import_glori()

########################################################################################################################
### GLORI p-val versus stoichiometry ###################################################################################
########################################################################################################################
# num_bins = 100
# range_Pvalue = [0, 0.002]
# range_Ratio = [0, 100]
# bin_counts, bin_x, bin_y = np.histogram2d(np.log10(df_glori['Pvalue']+1), df_glori['Ratio'], bins=num_bins, range=[range_Pvalue, range_Ratio])
# fig_glori_pval_ratio = plt.figure(figsize=(5, 5))
# plt.imshow(bin_counts, origin='lower')

thresh_p = 0.001
num_bins = 50
bin_range = [0, 1]
df_glori_low_pval = df_glori[df_glori['Pvalue']<thresh_p]
df_glori_high_pval = df_glori[df_glori['Pvalue']>=thresh_p]

plt.figure(figsize=(5*cm, 5*cm))
plt.subplot(2, 1, 1)
plt.hist(df_glori_low_pval['Ratio'], range=bin_range, bins=num_bins, label=f'Pvalue<{thresh_p}\n{len(df_glori_low_pval)} sites')
plt.xticks(np.arange(0, 1.01, 0.25))
plt.ylabel('Site Counts')
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.hist(df_glori_high_pval['Ratio'], range=bin_range, bins=num_bins, label=f'Pvalue$\geq${thresh_p}\n{len(df_glori_high_pval)} sites')
plt.xticks(np.arange(0, 1.01, 0.25))
plt.xlabel('$S_{GLORI}$')
plt.ylabel('Site Count')
plt.legend(loc='upper right')

plt.savefig(os.path.join(img_out, f'hist_GLORI_Pvalue_S.{FMT}'), **fig_kwargs)

# print(len(df_glori_high_pval) / len(df_glori))

########################################################################################################################
### scatter plot with GLORI high p-val #################################################################################
########################################################################################################################
df_merged = pd.merge(df_glori_high_pval, df_mAFiA, on=['chrom', 'chromStart', 'ref5mer'], suffixes=['_glori', '_mafia'])
df_merged_thresh = df_merged[df_merged['coverage']>=50]

# num_bins = 20
# counts, bins_x, bins_y = np.histogram2d(df_merged_thresh['modRatio_mafia'], df_merged_thresh['modRatio_glori'], bins=num_bins, range=[[0, 100], [0, 100]])

plt.figure(figsize=(4*cm, 4*cm))
plt.scatter(df_merged_thresh['modRatio_glori'], df_merged_thresh['modRatio_mafia'], s=1)
# plt.imshow(counts, origin='lower')
plt.xlim([0, 100])
plt.ylim([0, 100])
plt.xlabel('$S_{GLORI}$')
plt.ylabel('$S_{mAFiA}$')
plt.title(f'GLORI Pvalue $\geq$ {thresh_p}')

plt.savefig(os.path.join(img_out, f'scatter_mAFiA_GLORI_high_Pvalue.{FMT}'), **fig_kwargs)

########################################################################################################################
### histogram stoichiometry vs. coverage ###############################################################################
########################################################################################################################
coverage = df_mAFiA['coverage'].values
stoichio = df_mAFiA['modRatio'].values

num_bins = 20
coverage_range = [0, 500]
stoichio_range = [0, 100]

fig_cov_stoichio = plt.figure(figsize=(6, 6))
gs = fig_cov_stoichio.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
ax = fig_cov_stoichio.add_subplot(gs[1, 0])
ax_histx = fig_cov_stoichio.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig_cov_stoichio.add_subplot(gs[1, 1], sharey=ax)
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histy.tick_params(axis="y", labelleft=False)

bin_counts, bin_x, bin_y = np.histogram2d(coverage, stoichio, bins=num_bins, range=[coverage_range, stoichio_range])
ax.imshow(np.log10(bin_counts+1), origin='lower')
ax.set_xticks(np.arange(len(bin_x))[::2]-0.5, np.int32(bin_x)[::2])
ax.set_yticks(np.arange(len(bin_y))[::2]-0.5, np.int32(bin_y)[::2])
ax.set_xlabel('Coverage')
ax.set_ylabel('Stoichiometry')

x_bin_counts, _ = np.histogram(coverage, bins=bin_x)
ax_histx.bar(np.arange(len(x_bin_counts)), x_bin_counts)
ax_histx.set_xticks(np.arange(len(bin_x))[::2]-0.5, np.int32(bin_x)[::2])
ax_histx.set_ylabel('Num. Sites')

y_bin_counts, _ = np.histogram(stoichio, bins=bin_y)
ax_histy.barh(np.arange(len(y_bin_counts)), y_bin_counts)
ax_histy.set_yticks(np.arange(len(bin_y))[::2]-0.5, np.int32(bin_y)[::2])
ax_histy.set_xlabel('Num. Sites')

fig_cov_stoichio.savefig(os.path.join(img_out, f'coverage_stoichiometry.{FMT}'), **fig_kwargs)

fig_hist_stoichio, ax = plt.subplots(nrows=1, ncols=1, figsize=(5*cm, 5*cm))
ax.hist(df_mAFiA['modRatio'].values, range=[0, 100], bins=20)
ax.set_xlim(0, 100)
ax.set_xlabel('$S_{mAFiA}$')
ax.set_ylabel('Num. Sites')
fig_hist_stoichio.savefig(os.path.join(img_out, f'hist_stoichio.{FMT}'), **fig_kwargs)


########################################################################################################################
### correlation vs. coverage ###########################################################################################
########################################################################################################################
coverage_corr = []
for this_thresh_cov in range(10, 100):
    this_df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage'] >= this_thresh_cov]
    this_df_merged = pd.merge(this_df_mAFiA_thresh, df_glori, on=['chrom', 'chromStart'], suffixes=('_mafia', '_glori'))
    this_corr = np.corrcoef(this_df_merged[['modRatio_mafia', 'modRatio_glori']].values.T)[0, 1]
    coverage_corr.append((this_thresh_cov, this_corr))
coverage_corr = np.vstack(coverage_corr).T

fig_cov_corr, ax = plt.subplots(nrows=1, ncols=1, figsize=(5*cm, 5*cm))
ax.plot(coverage_corr[0], coverage_corr[1])
ax.set_xlabel('Min. Coverage')
ax.set_ylabel('$C_{(GLORI,mAFiA)}$')
fig_cov_corr.savefig(os.path.join(img_out, f'cov_corr.{FMT}'), **fig_kwargs)

########################################################################################################################
### correlation vs. stoichiometry ######################################################################################
########################################################################################################################
thresh_coverage = 50
df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage'] >= thresh_coverage]
df_merged = pd.merge(df_mAFiA_thresh, df_glori, on=['chrom', 'chromStart'], suffixes=('_mafia', '_glori'))

bin_width = 10
stoichio_lbins = np.arange(10, 100, bin_width)
bin_centers = (stoichio_lbins + stoichio_lbins + bin_width) / 2

stoichio_rms = []
for this_lbin in stoichio_lbins:
    sub_df_merged = df_merged[(df_merged['modRatio_glori']>=this_lbin) * (df_merged['modRatio_glori']<(this_lbin+bin_width))]
    # mafia_normed = np.float64(sub_df_merged['modRatio_mafia'].values)
    # mafia_normed -= mafia_normed.mean()
    # mafia_normed /= mafia_normed.std()
    # glori_normed = np.float64(sub_df_merged['modRatio_glori'].values)
    # glori_normed -= glori_normed.mean()
    # glori_normed /= glori_normed.std()
    # stoichio_corr.append(np.corrcoef(mafia_normed, glori_normed)[0, 1])
    stoichio_rms.append(np.sqrt(((sub_df_merged['modRatio_mafia'] - sub_df_merged['modRatio_glori'])**2).mean()))

stoichio_rms_norm = np.array(stoichio_rms) / stoichio_lbins

xticks = list(stoichio_lbins)
xticks.append(stoichio_lbins[-1]+bin_width)

fig_stoichio_rms, ax = plt.subplots(nrows=1, ncols=1, figsize=(5*cm, 5*cm))
ax.plot(bin_centers, stoichio_rms_norm)
ax.set_xticks(xticks)
ax.set_xlim([10, 100])
ax.set_xlabel('$S_{GLORI}$')
ax.set_ylabel('${\Delta}S_{(GLORI,mAFiA)}$ / $S_{GLORI}$')
fig_stoichio_rms.savefig(os.path.join(img_out, f'stoichio_RMS.{FMT}'), **fig_kwargs)


########################################################################################################################
### correlation ########################################################################################################
########################################################################################################################
thresh_coverage = 50
df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage'] >= thresh_coverage]
df_merged = pd.merge(df_mAFiA_thresh, df_glori, on=['chrom', 'chromStart', 'ref5mer'], suffixes=('_mafia', '_glori'))

df_merged = df_merged[[
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
    'coverage',
    'modRatio_mafia',
    'modRatio_glori',
    'Gene',
    'Transcript',
    'P_adjust'
]]
# df_merged.to_csv(os.path.join(img_out, 'df_merged_mAFiA_GLORI.tsv'), sep='\t', index=False)
df_merged.to_csv(os.path.join('/home/adrian/NCOMMS_revision/source_data/HEK293', 'source_data_Figure_2d.tsv'), sep='\t', index=False)


# df_merged_sel = df_merged[df_merged['P_adjust'] < 1E-10]
df_merged_sel = df_merged
motif_counts = Counter(df_merged_sel['ref5mer']).most_common()
motifs = [pair[0] for pair in motif_counts]

total_num_sites = df_merged_sel.shape[0]
total_corr = np.corrcoef(df_merged_sel['modRatio_glori'], df_merged_sel['modRatio_mafia'])[0, 1]

ordered_motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']

df_merged_6motifs = df_merged_sel[df_merged_sel['ref5mer'].isin(ordered_motifs)]
num_sites_6motifs = df_merged_6motifs.shape[0]
corr_6motifs = np.corrcoef(df_merged_6motifs['modRatio_glori'], df_merged_6motifs['modRatio_mafia'])[0, 1]

with open(os.path.join(img_out, 'sites_correlation_6motifs.txt'), 'w') as f_corr:
    f_corr.writelines(f'{num_sites_6motifs} sites, corr. {corr_6motifs:.2f}')

print(f'{num_sites_6motifs} sites, corr. {corr_6motifs:.2f}')
# num_motifs = 6
# num_rows = 2
# num_cols = 3

def scatter_plot_mafia_vs_glori(num_motifs, num_rows, num_cols, plot_name, figsize):
    ordered_motifs = [k for k, v in Counter(df_merged['ref5mer']).most_common()][:num_motifs]
    fig_corr, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=figsize)
    for ind, motif in enumerate(ordered_motifs):
        ax = axs[ind // num_cols, ind % num_cols]
        sub_df = df_merged_sel[df_merged_sel['ref5mer'] == motif]
        motif_corr = np.corrcoef(sub_df['modRatio_glori'], sub_df['modRatio_mafia'])[0, 1]
        # plt.subplot(num_rows, num_cols, ind + 1)
        # label = f"{motif.replace('T', 'U')}\n{corr:.2f}"
        label = f"{motif.replace('T', 'U')}\n({motif_corr:.2f})"
        if len(sub_df):
            ax.scatter(sub_df['modRatio_glori'], sub_df['modRatio_mafia'], s=0.1)
        # ax.set_title(f'{motif}, {corr:.2f}', x=0.26, y=1.0, pad=-15, backgroundcolor='black', color='white')
        # ax.set_title(f"{motif.replace('T', 'U')}", y=0.85)
        ax.set_xlim([-5, 105])
        ax.set_ylim([-5, 105])
        # ax.legend(loc='upper left', handlelength=0.1)
        ax.text(0, 80, label, fontsize=5)

        ax.set_xticks(range(0, 101, 25))
        ax.set_yticks(range(0, 101, 25))

        if ind%num_cols == 0:
            ax.set_yticks(range(0, 101, 25))
        else:
            ax.set_yticks([])
        if ind>=(num_rows-1)*num_cols:
            ax.set_xticks(range(0, 101, 25))
        else:
            ax.set_xticks([])

        # if ind%num_cols == 0:
        #     ax.set_ylabel('mAFiA')
        # if ind>=(num_rows-1)*num_cols:
        #     ax.set_xlabel('GLORI')
    # fig_corr.suptitle(f'HEK293 WT\nTotal {total_num_sites} sites\nCorr. {total_corr:.2f}')
    fig_corr.savefig(os.path.join(img_out, plot_name), **fig_kwargs)

scatter_plot_mafia_vs_glori(num_motifs=6, num_rows=2, num_cols=3, plot_name=f'corr_mAFiA_GLORI_DRACH_6motifs.{FMT}', figsize=(7.5*cm, 5*cm))
# scatter_plot_mafia_vs_glori(num_motifs=12, num_rows=3, num_cols=4, plot_name=f'corr_mAFiA_GLORI_DRACH_12motifs.{FMT}', figsize=(10*cm, 8*cm))
# scatter_plot_mafia_vs_glori(num_motifs=18, num_rows=3, num_cols=6, plot_name=f'corr_mAFiA_GLORI_DRACH_18motifs.{FMT}', figsize=(15*cm, 8*cm))