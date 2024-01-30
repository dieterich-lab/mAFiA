import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

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
# FMT = 'pdf'
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

img_out = '/home/adrian/img_out/Mainz_B01'
os.makedirs(img_out, exist_ok=True)

datasets = [
    'JK_HEK293_DMSO_merged',
    'JK_HEK293_STM2457_merged'
]
sel_chrom = '6'
sel_chromStart = 159724693
sel_chromEnd = 159758319

# datasets = [
#     'mESC_WT_DMSO_merged',
#     'mESC_Mettl3_KO_merged'
# ]
# sel_chrom = '17'
# sel_chromStart = 12964799
# sel_chromEnd = 12994543

suffixes = [ds.lstrip('JK_').rstrip('_merged') for ds in datasets]
display_suffixes = [suffix.replace('_', '-') for suffix in suffixes]

# bed_dmso = '/home/adrian/Data/TRR319_RMaP/Project_B01/Adrian/JK_HEK293_DMSO_merged/mAFiA/chrALL.mAFiA.sites.bed'
# df_dmso = pd.read_csv(bed_dmso, sep='\t', dtype={'chrom': str})
#
# bed_stm = '/home/adrian/Data/TRR319_RMaP/Project_B01/Adrian/JK_HEK293_STM2457_merged/mAFiA/chrALL.mAFiA.sites.bed'
# df_stm = pd.read_csv(bed_stm, sep='\t', dtype={'chrom': str})

dfs = [
    pd.read_csv(f'/home/adrian/Data/TRR319_RMaP/Project_B01/Adrian/{ds}/mAFiA/chrALL.mAFiA.sites.bed', sep='\t', dtype={'chrom': str})
    for ds in datasets
]

df_merged = pd.merge(*dfs,
                     on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                     # on=['chrom', 'chromStart', 'chromEnd', 'name', 'score'],
                     suffixes=[f'_{suf}' for suf in suffixes])

# plt.figure(figsize=(8, 8))
# plt.plot(df_merged['modRatio_dmso'], df_merged['modRatio_stm'], '.')
# plt.xlim([-1, 101])
# plt.ylim([-1, 101])

num_bins = 50
# vmax = 3
vmax = 30
ticks = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)
counts, bin_x, bin_y = np.histogram2d(
    df_merged[f'modRatio_{suffixes[1]}'], df_merged[f'modRatio_{suffixes[0]}'],
    bins=[num_bins, num_bins], range=[[0, 100], [0, 100]],
)
counts_log1p = np.log10(counts + 1)

fig_hist2d = plt.figure(figsize=(5*cm, 4.5*cm))
ax_hist2d = fig_hist2d.add_subplot()
im = ax_hist2d.imshow(counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
# im = ax_hist2d.imshow(counts_log1p, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
# ax_hist2d.axhline(y=np.where(bin_y==50)[0][0]-0.5, c='r', linestyle='--')
ax_hist2d.plot([0, num_bins-1], [0, num_bins-1], c='r', linestyle='--', alpha=0.5)
ax_hist2d.set_xticks(np.linspace(0, num_bins, 5)-0.5, ticks)
ax_hist2d.set_yticks(np.linspace(0, num_bins, 5)-0.5, ticks)
# cbar = fig_hist2d.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal', location='top')
cbar = fig_hist2d.colorbar(im, fraction=0.046, pad=0.04, orientation='vertical', location='right')
cbar.set_ticks(np.linspace(0, vmax, 3))
cbar.set_label('Site Density', rotation=-90, labelpad=10)
ax_hist2d.set_xlabel(suffixes[0])
ax_hist2d.set_ylabel(suffixes[1])
fig_hist2d.savefig(os.path.join(img_out, f'hist2d_{suffixes[0]}_{suffixes[1]}.{FMT}'), **fig_kwargs)

### WTAP gene ###
df_gene = df_merged[(df_merged['chrom']==sel_chrom) * (df_merged['chromStart']>=sel_chromStart) * (df_merged['chromEnd']<=sel_chromEnd)]
df_gene.to_csv(os.path.join(img_out, f'WTAP_{suffixes[0]}_{suffixes[1]}.tsv'), sep='\t', index=False)