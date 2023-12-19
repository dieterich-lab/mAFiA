import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
# import seaborn as sn
import numpy as np

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
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
fig_kwargs = dict(format=FMT, dpi=1200)
#######################################################################

df_all_aucs = pd.DataFrame([
    [0.95,	0.94,	0.94],
    [0.89,	0.96,	0.92],
    [0.94,	0.96,	0.96]],
    index=['RL', 'SL', 'RL+SL'],
    columns=['RL', 'SL', 'RL+SL'])

margin = 0.25
min_z = 0.85
xx = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2]) + margin
yy = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2]) + margin
zz = np.zeros(9) + min_z
dx = np.ones(9) - margin
dy = np.ones(9) - margin
dz = df_all_aucs.values.flatten() - min_z

fig_cross_validation = plt.figure(figsize=(12*cm, 10*cm))
ax = fig_cross_validation.add_subplot(1, 1, 1, projection='3d')
# sn.heatmap(ax=ax, data=df_all_aucs, annot=True, fmt='.2f', cbar=False)
ax.bar3d(xx, yy, zz, dx, dy, dz)
for (this_x, this_y, this_z) in zip(xx + 0.3, yy + 0.15, dz+min_z):
    ax.text(this_x, this_y, this_z, this_z, c='w')
ax.set_xticks([0.5, 1.5, 2.5], ['RL', 'SL', 'RL+SL'])
ax.xaxis.set_rotate_label(False)
ax.set_xlabel('Train', labelpad=0)
ax.set_yticks([0.5, 1.5, 2.5], ['RL', 'SL', 'RL+SL'])
ax.yaxis.set_rotate_label(False)
ax.set_ylabel('Validate', labelpad=0)
ax.set_zlim([min_z, 1.00])
ax.set_zticks(np.arange(min_z, 1.01, 0.05))
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('Average AUC', labelpad=0)
fig_cross_validation.savefig(os.path.join(img_out, f'all_validation_AUCs.{FMT}'), **fig_kwargs)
plt.close('all')