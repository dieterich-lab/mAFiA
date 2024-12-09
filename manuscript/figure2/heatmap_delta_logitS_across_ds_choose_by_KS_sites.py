import os
import pandas as pd
from functools import reduce
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
#######################################################################
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy
from collections import Counter

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
ks_file = os.path.join(res_dir, 'KS_test/annotated_diff_sites_TAC_HFpEF.tsv')
in_dir = os.path.join(res_dir, 'transcript_logit')
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

# ds = ['TAC', 'HFpEF', 'Diet']
ds = ['TAC', 'HFpEF']

ds_display = {
    'TAC': 'HFrEF',
    'HFpEF': 'HFpEF',
    # 'Diet': 'WD /\nCD'
}
ds_color = {
    'TAC': 'g',
    'HFpEF': 'r',
    'Diet': 'b'
}
mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

# thresh_pval = 0.05
num_gene_mods = 50

df_ks_sites = pd.read_csv(ks_file, sep='\t')

dfs = {}
for this_ds in ds:
    this_df = pd.read_csv(os.path.join(in_dir, f'delta_logitS_{this_ds}.tsv'), sep='\t')
    # this_df = this_df[this_df.pval < thresh_pval]
    this_df.rename(columns={
        'num_sites': f'num_sites_{this_ds}',
        'delta_logit': f'delta_logit_{this_ds}',
        'pval': f'pval_{this_ds}'
    }, inplace=True)
    dfs[this_ds] = this_df

df_merged = reduce(lambda left, right: pd.merge(left, right, on=['gene', 'mod']), list(dfs.values()))

########################################################################################################################

vmax = 0.5

fig = plt.figure(figsize=(8*cm, 15*cm))
for mod_ind, this_mod in enumerate(mods):
    df_ks_sites_mod = df_ks_sites[df_ks_sites['name'] == this_mod]
    sel_genes = [this_gene for (this_gene, this_count) in Counter(df_ks_sites_mod['gene']).most_common(num_gene_mods)]

    df_sel = df_merged[(df_merged['mod'] == this_mod) * df_merged['gene'].isin(sel_genes)]
    df_sel.drop_duplicates(inplace=True)
    df_sel.sort_values('delta_logit_TAC', ascending=False, inplace=True)

    sorted_gene = df_sel['gene'].values
    sorted_mat_delta_logit = df_sel[[f'delta_logit_{this_ds}' for this_ds in ds]].values

    plt.subplot(1, 2, mod_ind+1)
    im = plt.imshow(sorted_mat_delta_logit, aspect=0.25, vmin=-vmax, vmax=vmax, cmap='seismic')
    plt.xticks(np.arange(len(ds)), list(ds_display.values()))
    plt.yticks(np.arange(num_gene_mods), sorted_gene)
    plt.gca().xaxis.tick_top()

    plt.title(f'${{{dict_mod_display[this_mod]}}}$')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.2, 0.05, 0.6])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_ticks(np.linspace(-vmax, vmax, 5))
plt.savefig(os.path.join(img_out, f'heatmap_delta_logit_S.{FMT}'), **fig_kwargs)
