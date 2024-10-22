import os
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
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

protein_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'
IR_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/IRFinder/gene_log2fc_IRratio_pval_TAC.tsv'
polyA_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA/gene_polyA_log2fc_pval_TAC.tsv'

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'

thresh_pval = 0.05
thresh_num_IR = 10
thresh_num_reads = 20
df_protein = pd.read_csv(protein_file, sep='\t')
df_protein = df_protein[df_protein['tp'] == 'Combined']
df_protein = df_protein[df_protein['p.adj'] < thresh_pval]
df_protein.rename(columns={'SYMBOL': 'gene', 'logFC': 'log2fc'}, inplace=True)
df_IR = pd.read_csv(IR_file, sep='\t', dtype={'chrom': str})
df_IR = df_IR[df_IR['num_IR_junctions'] >= thresh_num_IR]
# df_IR = df_IR[df_IR['pval'] < thresh_pval]
df_polyA = pd.read_csv(polyA_file, sep='\t')
df_polyA = df_polyA[
    (df_polyA['pval'] < thresh_pval)
    * (df_polyA['num_reads_SHAM'] >= thresh_num_reads)
    * (df_polyA['num_reads_TAC'] >= thresh_num_reads)
    ]

def get_binned_vec_xy(df_x, df_y, bin_max, bin_width):
    common_genes = list(set(df_x['gene']).intersection(set(df_y['gene'])))
    log2fc_x = df_x.set_index('gene').loc[common_genes]['log2fc'].values
    log2fc_y = df_y.set_index('gene').loc[common_genes]['log2fc'].values

    # bin_max = 1.0
    # bin_width = 0.5
    num_bins = int(2 * bin_max / bin_width)
    bin_edges = np.linspace(-bin_max, bin_max, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    binned_y = []
    for bin_i in range(len(bin_centers)):
        bin_start = bin_edges[bin_i]
        bin_end = bin_edges[bin_i + 1]
        binned_y.append(log2fc_y[(log2fc_x >= bin_start) * (log2fc_x < bin_end)])
    mean_y = np.array([np.mean(ll) for ll in binned_y])
    # err_y = np.array([np.std(ll) for ll in binned_y])
    err_y = np.array([0.25 * (np.quantile(ll, 0.75) - np.quantile(ll, 0.25)) for ll in binned_y])

    return bin_centers, mean_y, err_y, bin_edges

# common_genes = list(set(df_protein['SYMBOL']).intersection(set(df_polyA['gene'])))
# log2fc_protein = df_protein.set_index('SYMBOL').loc[common_genes]['logFC'].values
# log2fc_polyA = df_polyA.set_index('gene').loc[common_genes]['log2fc'].values

# common_genes = list(set(df_protein['SYMBOL']).intersection(set(df_IR['gene'])))
# log2fc_IR = df_IR.set_index('gene').loc[common_genes]['log2fc'].values

# vec_x = log2fc_polyA
# vec_x = log2fc_IR
# vec_y = log2fc_protein

x_max = 1
x_bin_width = 0.5
vec_IR_protein = get_binned_vec_xy(df_IR, df_protein, bin_max=x_max, bin_width=x_bin_width)
vec_polyA_protein = get_binned_vec_xy(df_polyA, df_protein, bin_max=x_max, bin_width=x_bin_width)

y_bin_width = 0.1

fig = plt.figure(figsize=(4*cm, 4*cm))
# plt.boxplot(binned_y, positions=bin_centers, widths=bin_width*0.5, showfliers=False)

plt.subplot(2, 1, 1)
# ymax = np.round(vec_IR_protein[1].max(), 2)
# ymin = np.round(vec_IR_protein[1].min(), 2)
# yticks = np.round(np.linspace(ymin, ymax, 5), 2)
ymax = np.round((vec_IR_protein[1] + vec_IR_protein[2]).max() / y_bin_width) * y_bin_width
ymin = np.round((vec_IR_protein[1] - vec_IR_protein[2]).min() / y_bin_width) * y_bin_width
yticks = np.round(np.arange(ymin, ymax+0.01, y_bin_width), 2)
plt.plot(vec_IR_protein[0], vec_IR_protein[1], '-', c='b', label='IR ratio')
# plt.plot(vec_IR_protein[0], vec_IR_protein[1], '.', markersize=3, c='b')
plt.errorbar(vec_IR_protein[0], vec_IR_protein[1], marker='.', color='b', yerr=vec_IR_protein[2], capsize=3.0)
plt.xlim([-x_max, x_max])
plt.xticks(vec_IR_protein[3], vec_IR_protein[3])
plt.tick_params('x', top=True, labeltop=True, bottom=False, labelbottom=False, labelcolor='b')
plt.yticks(yticks, yticks)

plt.subplot(2, 1, 2)
# ymax = np.round(vec_polyA_protein[1].max(), 1)
# ymin = np.round(vec_polyA_protein[1].min(), 1)
# yticks = np.round(np.linspace(ymin, ymax, 5), 2)
ymax = np.round((vec_polyA_protein[1] + vec_polyA_protein[2]).max() / y_bin_width) * y_bin_width
ymin = np.round((vec_polyA_protein[1] - vec_polyA_protein[2]).min() / y_bin_width) * y_bin_width
yticks = np.round(np.arange(ymin, ymax+0.01, y_bin_width), 2)
plt.plot(vec_polyA_protein[0], vec_polyA_protein[1], '-', c='m', label='polyA len')
# plt.plot(vec_polyA_protein[0], vec_polyA_protein[1], '.', markersize=3, c='m')
plt.errorbar(vec_polyA_protein[0], vec_polyA_protein[1], marker='.', color='m', yerr=vec_polyA_protein[2], capsize=3.0)
plt.xticks(vec_polyA_protein[3], vec_polyA_protein[3])
plt.tick_params('x', labelcolor='m')
plt.yticks(yticks, yticks)

fig.savefig(os.path.join(img_out, f'mean_log2fc_protein_versus_IR_polyA.{FMT}'), **fig_kwargs)