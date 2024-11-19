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

thresh_pval = 0.05
thresh_log_pval = -np.log10(thresh_pval)
thresh_num_sites = 1
thresh_coverage = 1

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
polyA_file = os.path.join(base_dir, f'polyA/gene_polyA_log2fc_pval_{ds}.tsv')

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'
os.makedirs(img_out, exist_ok=True)

df_polyA = pd.read_csv(polyA_file, sep='\t')
# df_polyA = df_polyA[
#     (df_polyA['pval'] < thresh_pval)
#     # * (df_polyA[f"num_reads_{conditions[0].rstrip('_merged')}"] >= thresh_coverage)
#     # * (df_polyA[f"num_reads_{conditions[1].rstrip('_merged')}"] >= thresh_coverage)
#     ]
df_polyA.sort_values('pval', inplace=True)

vec_gene, vec_log2fc, vec_pval = df_polyA[['gene', 'log2fc', 'pval']].values.T
vec_log2fc = np.float64(vec_log2fc)
vec_log_pval = -np.log10(np.float64(vec_pval))
mask_up = (vec_log_pval >= thresh_log_pval) * (vec_log2fc >= 0)
mask_down = (vec_log_pval >= thresh_log_pval) * (vec_log2fc < 0)

xmax = 5
ymax = 6
xticks = np.linspace(-xmax, xmax, 5)
pt_size = 1
margin = 0.95

plt.figure(figsize=(4*cm, 4*cm))
plt.scatter(vec_log2fc, vec_log_pval, s=pt_size, c='gray', alpha=0.5, rasterized=True, edgecolors='none')
plt.scatter(vec_log2fc[mask_up], vec_log_pval[mask_up], s=pt_size, c='red', rasterized=True, edgecolors='none')
plt.scatter(vec_log2fc[mask_down], vec_log_pval[mask_down], s=pt_size, c='blue', rasterized=True, edgecolors='none')
plt.xlim([-xmax, xmax])
plt.ylim([0, ymax])
plt.xticks(xticks)
# plt.xlabel(f'log2fc transcript ${{{dict_mod_display[this_mod]}}}$')
# plt.ylabel('-log10(pval)')

plt.axhline(y=thresh_log_pval, c='gray', ls='--')

display_genes = 10

reg_genes = {}
reg_genes['up'] = vec_gene[mask_up]
reg_genes['down'] = vec_gene[mask_down]
reg_genes['up'] = [this_gene for this_gene in vec_gene[mask_up] if this_gene[:2] not in ['Gm', 'mt']]
reg_genes['down'] = [this_gene for this_gene in vec_gene[mask_down] if this_gene[:2] not in ['Gm', 'mt']]
plt.text(-xmax*margin, ymax*margin, '\n'.join(reg_genes['down'][:display_genes]), c='blue', ha='left', va='top')
plt.text(xmax*margin, ymax*margin, '\n'.join(reg_genes['up'][:display_genes]), c='red', ha='right', va='top')
plt.text(xmax*margin, ymax+0.1, f"{len(reg_genes['up'])} up", ha='right', va='bottom', c='red')
plt.text(-xmax*margin, ymax+0.1, f"{len(reg_genes['down'])} down", ha='left', va='bottom', c='blue')

plt.savefig(os.path.join(img_out, f'volcano_plot_polyA_{ds}_pval{thresh_pval}.{FMT}'), **fig_kwargs)