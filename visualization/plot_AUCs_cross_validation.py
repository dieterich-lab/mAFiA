import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sn

img_out = '/home/adrian/NCOMMS_revision/images/VALIDATION'

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

df_all_aucs = pd.DataFrame([
    [0.95,	0.94,	0.94],
    [0.89,	0.96,	0.92],
    [0.94,	0.96,	0.96]],
    index=['RL', 'SL', 'RL+SL'],
    columns=['RL', 'SL', 'RL+SL'])

fig_cross_validation = plt.figure(figsize=(3*cm, 6*cm))
ax = fig_cross_validation.subplots(1, 1)
sn.heatmap(ax=ax, data=df_all_aucs, annot=True, fmt='.2f', cbar=False)
fig_cross_validation.savefig(os.path.join(img_out, f'all_validation_AUCs.{FMT}'), **fig_kwargs)