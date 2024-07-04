import os
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

chrom = 'X'
# chromEnd = 56792649
# chromEnd = 56791998
days = ['day1', 'day7', 'day21']
isoforms = ['ENSMUST00000023854', 'ENSMUST00000114772']
colors = ['red', 'purple']
iso_colors = {iso: color for iso, color in zip(isoforms, colors)}

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

res_directory = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'

plt.figure(figsize=(5, 5))
for this_isoform in isoforms:
    day_values = []
    for this_day in days:
        df = pd.read_csv(os.path.join(res_directory, f'TAC_{this_day}', 'bambu', f'Fhl1.{this_isoform}.bed'), sep='\t', dtype={'chrom': str})
        df_sel = df[(df['chrom']==chrom)*(df['chromEnd']==chromEnd)]
        day_values.append(df_sel['modRatio'].values[0])
        mod = df_sel['name'].values[0]
    plt.plot(day_values, '-o', label=this_isoform, c=iso_colors[this_isoform])
    # plt.plot(day_values, 'o', label=this_isoform, c=iso_colors[this_isoform])
plt.ylim([60, 100])
plt.ylabel(f'$S_{{{dict_mod_display[mod]}}}$', fontsize=12)
plt.xticks(range(len(days)), days)
plt.title(f'chr{chrom}: {chromEnd}', fontsize=15)
