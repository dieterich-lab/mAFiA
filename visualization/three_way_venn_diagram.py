import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

bed6_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
ds = ['TAC', 'HFpEF', 'Diet']
# conditions = ['TAC', 'HFpEF', 'WT_WD']
conditions = ['SHAM', 'ctrl', 'WT_CD']

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

img_out = f"/home/adrian/img_out/venn_diagram_{'_'.join(ds)}"
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


plt.figure(figsize=(10, 5))
# this_mod = 'm6A'
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    sites = {}
    for this_ds in ds:
        this_df_mod = dfs[this_ds][dfs[this_ds]['name'] == this_mod]
        sites[this_ds] = set([tuple(val) for val in this_df_mod[bed6_fields].values])
    venn3(sites.values(), conditions)
    plt.title(fr'${{{dict_mod_display[this_mod]}}}$')
plt.suptitle(rf'$S\geq{thresh_modRatio}$')
plt.savefig(os.path.join(img_out, f"{'_'.join(conditions)}.png"), bbox_inches='tight')