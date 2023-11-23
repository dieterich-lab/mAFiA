import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from collections import Counter

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/HEK293/WT'
chrs = [str(i) for i in range(2, 8)] + ['18']
img_out = '/home/adrian/NCOMMS_revision/images/HEK293'
os.makedirs(img_out, exist_ok=True)

glori_path = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv'
miclip_path = '/home/adrian/Data/DRACH/miCLIP_union_flat_exclude_Y_chromosome.bed'

def import_mAFiA(thresh_coverage=10):
    dfs = []
    for this_chr in chrs:
        dfs.append(pd.read_csv(os.path.join(source_data_dir, f'chr{this_chr}', 'mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t'))
    df_mAFiA = pd.concat(dfs, ignore_index=True)
    df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage']>=thresh_coverage]
    return df_mAFiA_thresh

def import_glori(thresh_pval=1E-3):
    df_glori = pd.read_csv(glori_path, dtype={'Chr': str})
    df_glori['chrom'] = [chr.lstrip('chr') for chr in df_glori['Chr']]
    df_glori['chromStart'] = df_glori['Sites'] - 1
    df_glori['modRatio'] = np.int32(np.round(df_glori['NormeRatio'] * 100.0))
    df_glori_thresh = df_glori[df_glori['P_adjust']<thresh_pval]

    return df_glori_thresh

def import_miclip(thresh_source_counts=3):
    df_miclip = pd.read_csv(miclip_path, sep='\t',
                            names=['chrom', 'chromStart', 'chromEnd', 'source', 'score', 'strand'],
                            dtype={'chrom': str})
    df_miclip['modRatio'] = 100
    df_miclip['source_counts'] = [len(this_source.split(',')) for this_source in df_miclip['source'].values]
    df_miclip_thresh = df_miclip[df_miclip['source_counts']>=thresh_source_counts]
    return df_miclip_thresh

df_mAFiA = import_mAFiA()
df_glori = import_glori()
df_miclip = import_miclip()


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

# fig_hist_stoichio, ax = plt.subplots(nrows=1, ncols=1, figsize=(5*cm, 5*cm))
# ax.hist(df_mAFiA['modRatio'].values, range=[0, 100], bins=100)
# ax.set_xlim(0, 100)
# ax.set_xlabel('Stoichiometry')
# ax.set_ylabel('Counts')
# fig_hist_stoichio.savefig(os.path.join(img_out, f'hist_stoichio.{FMT}'), **fig_kwargs)


########################################################################################################################
### venn diagram #######################################################################################################
########################################################################################################################
thresh_coverage = 50
thresh_stoichio = 80

def draw_venn_diagram(dict_dfs, thresh_stoichiometry=80):
    name_sites = {}
    for name, df in dict_dfs.items():
        df_sel = df[(df['modRatio']>=thresh_stoichiometry)*(df['chrom'].isin(chrs))]
        name_sites[name] = set([tuple(val) for val in df_sel[['chrom', 'chromStart']].values])

        # num_intersection = len(set(df1_sites).intersection(set(df2_sites)))
        # print(len(df1_sites), len(df2_sites), num_intersection)

    return venn3(name_sites.values(), name_sites.keys())

df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage']>=thresh_coverage]

venn_dict = {
    'mAFiA': df_mAFiA_thresh,
    'miClip': df_miclip,
    'GLORI': df_glori
}

fig_venn, ax = plt.subplots(nrows=1, ncols=1, figsize=(10*cm, 10*cm))
v = draw_venn_diagram(venn_dict, thresh_stoichiometry=thresh_stoichio)
fig_venn.savefig(os.path.join(img_out, f'venn_diagram.{FMT}'), **fig_kwargs)


########################################################################################################################
### correlation ########################################################################################################
########################################################################################################################
df_merged = pd.merge(df_mAFiA_thresh, df_glori, on=['chrom', 'chromStart'], suffixes=('_mafia', '_glori'))
# df_merged_sel = df_merged[df_merged['P_adjust'] < 1E-10]
df_merged_sel = df_merged
motif_counts = Counter(df_merged_sel['ref5mer']).most_common()
motifs = [pair[0] for pair in motif_counts]

total_num_sites = df_merged_sel.shape[0]
total_corr = np.corrcoef(df_merged_sel['modRatio_glori'], df_merged_sel['modRatio_mafia'])[0, 1]

ordered_motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']
num_rows = 3
num_cols = 2

fig_corr, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(4*cm, 6*cm))
for ind, motif in enumerate(ordered_motifs):
    ax = axs[ind // num_cols, ind % num_cols]
    sub_df = df_merged_sel[df_merged_sel['ref5mer'] == motif]
    corr = np.corrcoef(sub_df['modRatio_glori'], sub_df['modRatio_mafia'])[0, 1]
    # plt.subplot(num_rows, num_cols, ind + 1)
    label = f"{motif.replace('T', 'U')}\n{corr:.2f}"
    ax.scatter(sub_df['modRatio_glori'], sub_df['modRatio_mafia'], s=0.1)
    # ax.set_title(f'{motif}, {corr:.2f}', x=0.26, y=1.0, pad=-15, backgroundcolor='black', color='white')
    # ax.set_title(f"{motif.replace('T', 'U')}", y=0.85)
    ax.set_xlim([-5, 105])
    ax.set_ylim([-5, 105])
    # ax.legend(loc='upper left', handlelength=0.1)
    ax.text(0, 75, label)

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
fig_corr.savefig(os.path.join(img_out, f'corr_mAFiA_GLORI_DRACH.{FMT}'), **fig_kwargs)
