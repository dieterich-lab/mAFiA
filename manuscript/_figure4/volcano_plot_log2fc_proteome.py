import pandas as pd
import numpy as np
import os
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

data_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/protein_abundance'
proteome_file = os.path.join(data_dir, 'TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv')

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/_figure4'

df_proteome = pd.read_csv(proteome_file, sep='\t')
df_proteome = df_proteome[df_proteome['tp'] == 'Combined']

xmax = 2
ymax = 7
xticks = np.linspace(-xmax, xmax, 5)
pt_size = 1
thresh_log_pval = 4
thresh_log2fc = 0.5

reg_genes = {}

plt.figure(figsize=(4*cm, 6*cm))
vec_gene, vec_log2fc, vec_log_pval = df_proteome[['SYMBOL', 'logFC', 'LogPVal']].values.T

mask_up = (vec_log_pval >= thresh_log_pval) * (vec_log2fc > thresh_log2fc)
mask_down = (vec_log_pval >= thresh_log_pval) * (vec_log2fc < -thresh_log2fc)

plt.scatter(vec_log2fc, vec_log_pval, s=pt_size, c='gray', alpha=0.5, rasterized=True, edgecolors='none')
plt.scatter(vec_log2fc[mask_up], vec_log_pval[mask_up], s=pt_size, c='red', rasterized=True, edgecolors='none')
plt.scatter(vec_log2fc[mask_down], vec_log_pval[mask_down], s=pt_size, c='blue', rasterized=True, edgecolors='none')
plt.xlim([-xmax, xmax])
plt.ylim([0, ymax])
plt.xticks(xticks)
# plt.xlabel(f'log2fc transcript ${{{dict_mod_display[this_mod]}}}$')
# plt.ylabel('-log10(pval)')

plt.axhline(y=thresh_log_pval, c='gray', ls='--')

reg_genes['up'] = vec_gene[mask_up]
reg_genes['down'] = vec_gene[mask_down]
reg_genes['up'] = [this_gene for this_gene in vec_gene[mask_up] if this_gene[:2]!='Gm']
reg_genes['down'] = [this_gene for this_gene in vec_gene[mask_down] if this_gene[:2]!='Gm']
plt.text(-xmax*0.9, ymax*0.8, '\n'.join(reg_genes['down']), c='blue', ha='left', va='top')
plt.text(xmax*0.55, ymax*0.8, '\n'.join(reg_genes['up']), c='red', ha='left', va='top')

plt.savefig(os.path.join(img_out,
                         f'volcano_plot_proteome_logPval{thresh_log_pval}_log2fc{thresh_log2fc}.{FMT}'),
            **fig_kwargs)

with open(os.path.join(data_dir, 'up_down_reg_protein_TAC_SHAM_combined.txt'), 'w') as f_out:
    f_out.write('#up_regulated\n')
    for this_gene in reg_genes['up']:
        f_out.write(f'{this_gene}\n')
    f_out.write('#down_regulated\n')
    for this_gene in reg_genes['down']:
        f_out.write(f'{this_gene}\n')
