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

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

stoichiometry_file = os.path.join(
    '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_stoichiometry',
    'transcript_stoichiometry_TAC_merged_vs_SHAM_merged.tsv'
)

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'

df_stoichiometry = pd.read_csv(stoichiometry_file, sep='\t')

xmax = 4
ymax = 6
xticks = np.linspace(-xmax, xmax, 5)

pt_size = 1

thresh_log_pval = 4
thresh_log2fc = 0.25
reg_genes = {mod: {} for mod in mods}
plt.figure(figsize=(9*cm, 4*cm))
for mod_ind, this_mod in enumerate(mods):
    this_mod_df = df_stoichiometry[df_stoichiometry['mod'] == this_mod]
    vec_gene, vec_log2fc, vec_pval = this_mod_df[['gene', 'log2fc', 'pval']].values.T
    vec_log_pval = -np.log10(np.float64(vec_pval))

    mask_up = (vec_log_pval >= thresh_log_pval) * (vec_log2fc > thresh_log2fc)
    mask_down = (vec_log_pval >= thresh_log_pval) * (vec_log2fc < -thresh_log2fc)

    plt.subplot(1, 2, mod_ind+1)
    plt.scatter(vec_log2fc, vec_log_pval, s=pt_size, c='gray', alpha=0.5, rasterized=True, edgecolors='none')
    plt.scatter(vec_log2fc[mask_up], vec_log_pval[mask_up], s=pt_size, c='red', rasterized=True, edgecolors='none')
    plt.scatter(vec_log2fc[mask_down], vec_log_pval[mask_down], s=pt_size, c='blue', rasterized=True, edgecolors='none')
    plt.xlim([-xmax, xmax])
    plt.ylim([0, ymax])
    plt.xticks(xticks)
    # plt.xlabel(f'log2fc transcript ${{{dict_mod_display[this_mod]}}}$')
    # plt.ylabel('-log10(pval)')

    plt.axhline(y=thresh_log_pval, c='gray', ls='--')

    reg_genes[this_mod]['up'] = vec_gene[mask_up]
    reg_genes[this_mod]['down'] = vec_gene[mask_down]
    reg_genes[this_mod]['up'] = [this_gene for this_gene in vec_gene[mask_up] if this_gene[:2]!='Gm']
    reg_genes[this_mod]['down'] = [this_gene for this_gene in vec_gene[mask_down] if this_gene[:2]!='Gm']
    plt.text(-xmax+0.5, ymax-0.1, '\n'.join(reg_genes[this_mod]['down']), c='blue', ha='left', va='top')
    plt.text(xmax-2, ymax-0.1, '\n'.join(reg_genes[this_mod]['up']), c='red', ha='left', va='top')

plt.savefig(os.path.join(img_out,
                         f'volcano_plot_stoichiometry_logPval{thresh_log_pval}_log2fc{thresh_log2fc}.{FMT}'),
            **fig_kwargs)
