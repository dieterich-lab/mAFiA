import os
import pandas as pd
from matplotlib_venn import venn3
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

bed6_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
ds = ['TAC', 'Diet', 'HFpEF']
conditions = ['TAC', 'WT_WD', 'HFpEF']
# conditions = ['SHAM', 'WT_CD', 'ctrl']

ds_colors = {
    this_ds: this_color for (this_ds, this_color) in zip(ds, ['g', 'b', 'r'])
}

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'
os.makedirs(img_out, exist_ok=True)

thresh_confidence = 0.0
thresh_modRatio = 50.0

dfs = {}
for this_ds, this_cond in zip(ds, conditions):
    this_df = pd.read_csv(os.path.join(res_dir, this_ds, f'{this_cond}_merged', 'chrALL.mAFiA.sites.bed'),
                             sep='\t', dtype={'chrom': str})
    dfs[this_ds] = this_df[
        (this_df['confidence'] >= thresh_confidence)
        * (this_df['modRatio'] >= thresh_modRatio)
    ]


plt.figure(figsize=(12*cm, 6*cm))
# this_mod = 'm6A'
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    sites = {}
    for this_ds in ds:
        this_df_mod = dfs[this_ds][dfs[this_ds]['name'] == this_mod]
        sites[this_ds] = set([tuple(val) for val in this_df_mod[bed6_fields].values])
    venn3(sites.values(), ds, set_colors=ds_colors.values())
    plt.title(fr'${{{dict_mod_display[this_mod]}}}$')
# plt.suptitle(rf'$S\geq{thresh_modRatio}$')
plt.savefig(os.path.join(img_out, f'three_way_venn_diagram.{FMT}'), **fig_kwargs)