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

in_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/transcript_logit'
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

ds = ['TAC', 'HFpEF', 'Diet']
ds_display = {
    'TAC': 'TAC /\nSHAM',
    'HFpEF': 'HFpEF /\nCtrl',
    'Diet': 'WD /\nCD'
}
mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresh_pval = 0.05
num_gene_mods = 20

dfs = {}
for this_ds in ds:
    this_df = pd.read_csv(os.path.join(in_dir, f'delta_logitS_{this_ds}.tsv'), sep='\t')
    this_df = this_df[this_df.pval < thresh_pval]
    this_df.rename(columns={
        'num_sites': f'num_sites_{this_ds}',
        'delta_logit': f'delta_logit_{this_ds}',
        'pval': f'pval_{this_ds}'
    }, inplace=True)
    dfs[this_ds] = this_df

df_merged = reduce(lambda left, right: pd.merge(left, right, on=['gene', 'mod']), list(dfs.values()))

vec_gene_mod = [f'{gene}, ${dict_mod_display[mod]}$' for gene, mod in df_merged[['gene', 'mod']].values]
mat_delta_logit = df_merged[[f'delta_logit_{this_ds}' for this_ds in ds]].values

### select rows ###
# vec_variance = np.sum(np.abs(mat_delta_logit), axis=1)
# vec_variance = np.sum(mat_delta_logit**2, axis=1)
vec_variance = df_merged[[f'num_sites_{this_ds}' for this_ds in ds]].values.sum(axis=1)

sel_indices = np.argpartition(vec_variance, -num_gene_mods)[-num_gene_mods:]

for row in [vec_gene_mod[ind] for ind in sel_indices]:
    print(row)

sel_gene_mod = [vec_gene_mod[this_ind] for this_ind in sel_indices]
sel_mat_delta_logit = mat_delta_logit[sel_indices]

# sorted_indices = np.lexsort((
#         sel_mat_delta_logit[:, 0],
#         # np.round(sel_mat_delta_logit[:, 1]),
#         # np.round(sel_mat_delta_logit[:, 2])
#     ))
sorted_indices = np.argsort(sel_gene_mod)[::-1]
sorted_mat_delta_logit = sel_mat_delta_logit[sorted_indices]
sorted_gene_mod = np.array(sel_gene_mod)[sorted_indices]

vmax = 0.5
plt.figure(figsize=(4*cm, 8.5*cm))
plt.imshow(sorted_mat_delta_logit, aspect='auto', origin='lower', vmin=-vmax, vmax=vmax, cmap='seismic')
plt.xticks(np.arange(len(ds)), list(ds_display.values()))
plt.yticks(np.arange(num_gene_mods), sorted_gene_mod)
plt.gca().xaxis.tick_top()
cbar = plt.colorbar()
cbar.set_ticks(np.linspace(-vmax, vmax, 5))
plt.savefig(os.path.join(img_out, f'heatmap_delta_logit_S_all_ds.{FMT}'), **fig_kwargs)