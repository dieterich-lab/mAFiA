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

def get_logit(in_val):
    return np.log2(in_val / (1.0 - in_val))

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']
# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
stoichiometry_file = os.path.join(base_dir, 'transcript_logit', f'delta_logitS_{ds}.tsv')
polyA_file = os.path.join(base_dir, f'polyA/gene_polyA_log2fc_pval_{ds}.tsv')

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'
os.makedirs(img_out, exist_ok=True)

thresh_pval = 0.05
thresh_num_sites = 1
thresh_coverage = 1
log2fc_bin_boundary = 0.1

df_stoichiometry = pd.read_csv(stoichiometry_file, sep='\t')
df_stoichiometry = df_stoichiometry[
    (df_stoichiometry['pval'] < thresh_pval)
    * (df_stoichiometry[f'num_sites'] >= thresh_num_sites)
    ]
df_polyA = pd.read_csv(polyA_file, sep='\t')
df_polyA = df_polyA[
    (df_polyA['pval'] < thresh_pval)
    * (df_polyA[f"num_reads_{conditions[0].rstrip('_merged')}"] >= thresh_coverage)
    * (df_polyA[f"num_reads_{conditions[1].rstrip('_merged')}"] >= thresh_coverage)
    ]

common_genes = list(set(df_stoichiometry['gene']).intersection(set(df_polyA['gene'])))
vec_log2fc_mod = {}
for this_mod in mods:
    df_stoichiometry_this_mod = df_stoichiometry[df_stoichiometry['mod'] == this_mod]
    vec_log2fc_mod[this_mod] = np.array(
        # [df_stoichiometry_this_mod.loc[this_gene]['log2fc']
        [df_stoichiometry_this_mod[df_stoichiometry_this_mod['gene'] == this_gene]['delta_logit'].values[0]
         if this_gene in df_stoichiometry_this_mod['gene'].values else np.nan
         for this_gene in common_genes]
    )
vec_log2fc_polyA = np.array([df_polyA.set_index('gene').loc[this_gene]['log2fc'] for this_gene in common_genes])

# ylim = [-0.5, 0.0]
ylim = [-0.5, 0.5]
yticks = np.round(np.linspace(*ylim, 3), 2)
plt.figure(figsize=(9*cm, 4*cm))
for mod_ind, this_mod in enumerate(mods):
    # binned_y = [
    #     [x for x in vec_log2fc_mod[this_mod][(vec_log2fc_polyA < -log2fc_bin_boundary)] if ~np.isnan(x) and ~np.isinf(x)],
    #     [x for x in vec_log2fc_mod[this_mod][(vec_log2fc_polyA >= log2fc_bin_boundary)] if ~np.isnan(x) and ~np.isinf(x)],
    # ]
    binned_y = [
        [x for x in vec_log2fc_polyA[(vec_log2fc_mod[this_mod] < -log2fc_bin_boundary)] if ~np.isnan(x) and ~np.isinf(x)],
        [x for x in vec_log2fc_polyA[(vec_log2fc_mod[this_mod] >= log2fc_bin_boundary)] if ~np.isnan(x) and ~np.isinf(x)],
    ]
    num_genes = [len(this_bin) for this_bin in binned_y]
    plt.subplot(1, 2, mod_ind+1)
    plt.boxplot(binned_y,
                positions=[-0.5, 0.5],
                widths=0.25,
                whis=0.5,
                showfliers=False)
    # plt.ylim(ylim)
    # plt.xticks([-0.5, 0.5],
    #            [f'<-{log2fc_bin_boundary} ({num_genes[0]})', rf'$\geq${log2fc_bin_boundary} ({num_genes[1]})'])
    plt.xticks([-0.5, 0.5],
               [rf'$\downarrow$ ({num_genes[0]})', rf'$\uparrow$ ({num_genes[1]})'])
    plt.yticks(yticks)
    # plt.xlabel(r'$\langle$$\Delta$logit(S)$\rangle$')
    # if mod_ind == 0:
    #     plt.ylabel('log2fc polyA')
    # plt.title(rf'${{{dict_mod_display[this_mod]}}}$')
    # if mod_ind == 0:
    #     plt.yticks(yticks, yticks)
    # else:
    #     plt.yticks(yticks, [])
# plt.savefig(os.path.join(img_out, f'boxplot_log2fc_three_prime_utr_stoichiometry_vs_log2fc_polyA_{ds}.{FMT}'),
#             **fig_kwargs)
# plt.savefig(os.path.join(img_out, f'boxplot_delta_logit_S_vs_log2fc_polyA_{ds}.{FMT}'),
#             **fig_kwargs)
plt.savefig(os.path.join(img_out, f'boxplot_log2fc_polyA_vs_delta_logit_S_{ds}_pval{thresh_pval}.{FMT}'),
            **fig_kwargs)