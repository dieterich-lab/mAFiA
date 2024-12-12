import os
from glob import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

datasets = [f'Dataset{ii}' for ii in range(1, 7)]
data_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Test'

reference = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/test_reference.fasta'
with open(reference, 'r') as h_ref:
    for line in h_ref.readlines():
        if line[0]!='>':
            ref = line

dict_dfs = {ds: pd.read_csv(os.path.join(data_dir, ds, 'Gmorah.sites.bed'), sep='\t') for ds in datasets}
dfs = [df.rename(columns={'coverage': f'coverage_{ds}', 'modRatio': f'modRatio_{ds}'}) for ds, df in dict_dfs.items()]
df_merged = dfs[0].copy()
for ii in range(1, 6):
    df_merged = pd.merge(df_merged, dfs[ii], on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'])

array_modRatio = df_merged[[f'modRatio_Dataset{ii}' for ii in range(1, 7)]].values
positions = df_merged['chromStart'].values

plt.figure(figsize=(8, 10))
plt.imshow(array_modRatio)
plt.xticks(range(6), range(1, 7))
plt.yticks(range(len(positions)), positions)
plt.xlabel('Dataset')
plt.ylabel('Position')

std_modRatio = np.std(array_modRatio, axis=1)
sel_indices = np.argpartition(std_modRatio, -3)[-3:]
sel_indices.sort()
sel_positions = positions[sel_indices]
sel_modRatio = array_modRatio[sel_indices]

plt.figure(figsize=(5, 10))
for subplot_ind, (pos, ratio) in enumerate(zip(sel_positions, sel_modRatio)):
    plt.subplot(3, 1, subplot_ind+1)
    ref5mer = ref[pos-2:pos+3]
    plt.plot(ratio, '-o', label=f'Postion {pos}\n{ref5mer}')
    plt.ylim([20, 70])
    plt.legend(loc='upper right', fontsize=10)
    plt.xticks(range(6), range(1, 7))
    plt.ylabel('Mod. Ratio', fontsize=12)
plt.xlabel('Dataset', fontsize=12)
plt.savefig(os.path.join(img_out, 'variation_modRatio_by_position.png'), bbox_inches='tight')
