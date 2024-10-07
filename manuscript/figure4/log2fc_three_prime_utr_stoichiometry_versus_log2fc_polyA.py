import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import numpy as np
from tqdm import tqdm
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

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

ds = 'TAC'
base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
stoichiometry_file = os.path.join(base_dir, 'three_prime_utr_stoichiometry',
                                  'three_prime_utr_stoichiometry_TAC_merged_vs_SHAM_merged.tsv')
polyA_file = os.path.join(base_dir, 'polyA/gene_polyA_log2fc_pval_TAC.tsv')

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'
os.makedirs(img_out, exist_ok=True)

thresh_pval = 0.05
df_stoichiometry = pd.read_csv(stoichiometry_file, sep='\t')
df_stoichiometry = df_stoichiometry[df_stoichiometry['pval'] < thresh_pval]
df_polyA = pd.read_csv(polyA_file, sep='\t')
df_polyA = df_polyA[df_polyA['pval'] < thresh_pval]

common_genes = list(set(df_stoichiometry['gene']).intersection(set(df_polyA['gene'])))
vec_log2fc_mod = {}
for this_mod in mods:
    df_stoichiometry_this_mod = df_stoichiometry[df_stoichiometry['mod'] == this_mod].set_index('gene')
    vec_log2fc_mod[this_mod] = np.array(
        [df_stoichiometry_this_mod.loc[this_gene]['log2fc']
         if this_gene in df_stoichiometry_this_mod.index else np.nan
         for this_gene in common_genes]
    )
vec_log2fc_polyA = np.array([df_polyA.set_index('gene').loc[this_gene]['log2fc'] for this_gene in common_genes])

boundary = 0.1
ylim = [-0.2, 0.4]
yticks = np.round(np.linspace(*ylim, 4), 2)
plt.figure(figsize=(5*cm, 4*cm))
for mod_ind, this_mod in enumerate(mods):
    binned_y = [
        [x for x in vec_log2fc_mod[this_mod][(vec_log2fc_polyA < -boundary)] if ~np.isnan(x) and ~np.isinf(x)],
        [x for x in vec_log2fc_mod[this_mod][(vec_log2fc_polyA >= boundary)] if ~np.isnan(x) and ~np.isinf(x)],
    ]
    plt.subplot(1, 2, mod_ind+1)
    plt.boxplot(binned_y, positions=[-0.5, 0.5], widths=0.25, whis=0.5, showfliers=False)
    # plt.ylim(ylim)
    plt.xticks([-0.5, 0.5], [f'< -{boundary}', rf'$\geq$ {boundary}'])
    if mod_ind == 0:
        plt.yticks(yticks, yticks)
    else:
        plt.yticks(yticks, [])
plt.savefig(os.path.join(img_out, f'boxplot_log2fc_three_prime_utr_stoichiometry_vs_log2fc_polyA_{ds}.{FMT}'),
            **fig_kwargs)
