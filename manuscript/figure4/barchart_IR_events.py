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
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'

ds = ['TAC', 'HFpEF', 'Diet']
cond = {
    'TAC': ['SHAM_merged', 'TAC_merged'],
    'HFpEF': ['ctrl_merged', 'HFpEF_merged'],
    'Diet': ['WT_CD_merged', 'WT_WD_merged'],
}

ds_cond_IR_counts = {this_ds: [] for this_ds in ds}
for this_ds in ds:
    for this_cond in cond[this_ds]:
        filename = os.path.join(img_out, f'df_irf_clean_{this_ds}_{this_cond}.tsv')
        with open(filename, "rb") as f:
            ds_cond_IR_counts[this_ds].append(sum(1 for _ in f))

plt.figure(figsize=(4*cm, 4*cm))
for ds_ind, this_ds in enumerate(ds):
    plt.bar(ds_ind-0.15, ds_cond_IR_counts[this_ds][0], width=0.2, fc='b')
    plt.bar(ds_ind+0.15, ds_cond_IR_counts[this_ds][1], width=0.2, fc='r')
plt.xticks(range(len(ds)), ds)
plt.savefig(os.path.join(img_out, f'barchart_IR_events.{FMT}'), **fig_kwargs)
