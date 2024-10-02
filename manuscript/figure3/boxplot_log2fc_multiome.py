import os
import pandas as pd
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


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}


thresh_pval = 0.05

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']
cond_names = [this_cond.rstrip('_merged') for this_cond in conditions]

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'
os.makedirs(img_out, exist_ok=True)

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_TAC_vs_SHAM_benjamini_hochberg.tsv'
df_proteome_all = pd.read_csv(proteome_file, sep='\t')
df_proteome = df_proteome_all[
    (df_proteome_all['tp'] == 'Combined')
    * (df_proteome_all['p.adj'] < thresh_pval)
]
df_proteome.rename(columns={'logFC': 'log2fc_proteome'}, inplace=True)

epitranscriptome_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_stoichiometry/transcript_stoichiometry_TAC_merged_vs_SHAM_merged.tsv'
df_epitranscriptome_all = pd.read_csv(epitranscriptome_file, sep='\t')
df_epitranscriptome = df_epitranscriptome_all[
    df_epitranscriptome_all['pval'] < thresh_pval
]
# df_epitranscriptome = df_epitranscriptome_all
# df_epitranscriptome = df_epitranscriptome_all[
#     (df_epitranscriptome_all['num_nts_SHAM_merged'] >= 1000)
#     * (df_epitranscriptome_all['num_nts_TAC_merged'] >= 1000)
# ]
df_epitranscriptome.rename(columns={'log2fc': 'log2fc_epitranscriptome'}, inplace=True)

transcriptome_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_abundance/transcript_abundance_TAC_merged_vs_SHAM_merged.tsv'
df_transcriptome = pd.read_csv(transcriptome_file, sep='\t')
df_transcriptome.rename(columns={'log2fc': 'log2fc_transcriptome'}, inplace=True)

gene_transcriptome_epitranscriptome_proteome = {}
for this_gene in tqdm(df_proteome['SYMBOL']):
    # if (this_gene in df_transcriptome['gene'].values) and (this_gene in df_epitranscriptome['gene'].values):
    if (this_gene in df_transcriptome['gene'].values):
        sub_df_epitranscriptome = df_epitranscriptome[df_epitranscriptome['gene'] == this_gene]
        log2fc_mods = {}
        for this_mod in mods:
            if this_mod in sub_df_epitranscriptome['mod'].values:
                log2fc_mods[this_mod] = sub_df_epitranscriptome[
                    sub_df_epitranscriptome['mod'] == this_mod
                    ]['log2fc_epitranscriptome'].values[0]
            else:
                log2fc_mods[this_mod] = None
        gene_transcriptome_epitranscriptome_proteome[this_gene] = (
            df_transcriptome[df_transcriptome['gene'] == this_gene]['log2fc_transcriptome'].values[0],
            log2fc_mods['m6A'],
            log2fc_mods['psi'],
            df_proteome[df_proteome['SYMBOL'] == this_gene]['log2fc_proteome'].values[0]
        )

vec_abundance, vec_m6A, vec_psi, vec_protein = np.vstack(gene_transcriptome_epitranscriptome_proteome.values()).T

### boxplot ###
def get_1d_binned_values(in_vec_x, in_vec_y, xmax, num_bins):
    in_vec_x[in_vec_x == None] = -np.inf

    x_bin_edges = np.linspace(-xmax, xmax, num_bins + 1)
    x_bin_centers = 0.5 * (x_bin_edges[1:] + x_bin_edges[:-1])
    binned_y = []
    for bin_i in range(len(x_bin_edges) - 1):
        bin_start = x_bin_edges[bin_i]
        bin_stop = x_bin_edges[bin_i + 1]
        mask = (in_vec_x >= bin_start) * (in_vec_x < bin_stop)
        binned_y.append(in_vec_y[mask])
    medians_y = [np.mean(ll) for ll in binned_y]
    return x_bin_edges, x_bin_centers, binned_y, medians_y

def do_boxplot(vec_x, vec_y, xmax, ymax):
    x_bin_edges, x_bin_centers, binned_y, medians_y = get_1d_binned_values(vec_x, vec_y, xmax, 4)
    plt.boxplot(binned_y, positions=x_bin_centers, showfliers=False, whis=0.5, medianprops={'linewidth': 0})
    plt.plot(x_bin_centers, medians_y, 'r-')
    plt.plot(x_bin_centers, medians_y, 'r+', markersize=3)
    plt.xticks(x_bin_edges, x_bin_edges)
    plt.yticks(yticks)
    plt.xlim([-xmax, xmax])
    plt.ylim([-ymax, ymax])

ymax = 0.55
yticks = np.linspace(-0.5, 0.5, 5)

plt.figure(figsize=(13*cm, 4*cm))

plt.subplot(1, 3, 1)
do_boxplot(vec_abundance, vec_protein, 2, ymax)
# plt.xlabel('log2fc transcript abundance')
# plt.ylabel('log2fc protein abundance')

plt.subplot(1, 3, 2)
xmax = 1
do_boxplot(vec_m6A, vec_protein, 1, ymax)
# plt.xlabel(f"log2fc transcript ${{{dict_mod_display['m6A']}}}$")
# plt.ylabel('log2fc protein abundance')

plt.subplot(1, 3, 3)
do_boxplot(vec_psi, vec_protein, 1, ymax)
# plt.xlabel(f"log2fc transcript ${{{dict_mod_display['psi']}}}$")
# plt.ylabel('log2fc protein abundance')

plt.savefig(os.path.join(img_out, f'boxplot_log2fc_multiome.{FMT}'), **fig_kwargs)