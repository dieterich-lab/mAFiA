import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ds = 'Dataset5'
res_dir = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Test/{ds}'
img_out = '/home/adrian/img_out/Gmorah'
os.makedirs(img_out, exist_ok=True)

dfs = [pd.read_csv(os.path.join(res_dir, f'Gmorah_batch{i}/mAFiA.sites.bed'), sep='\t') for i in range(5)]

df_combined = dfs[0][['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer']]
coverage_combined = np.vstack([df['coverage'] for df in dfs]).sum(axis=0)
modRatio_combined = np.int32(np.round(np.vstack([df['coverage'] * df['modRatio'] for df in dfs]).sum(axis=0) / coverage_combined))

df_combined['coverage'] = coverage_combined
df_combined['modRatio'] = modRatio_combined
df_combined.to_csv(os.path.join(res_dir, 'Gmorah.sites.bed'), sep='\t', index=False)

### plot ###
loc_names = [f'{x}_{y}' for (x, y) in zip(df_combined['chromStart'], df_combined['ref5mer'])]

plt.figure(figsize=(16, 9))
plt.bar(range(len(loc_names)), df_combined['modRatio'])
plt.xticks(range(len(loc_names)), loc_names, rotation=90)
plt.xlim([-1, len(loc_names)])
plt.ylim([50, 100])
plt.ylabel('Mod. Ratio', fontsize=20)
plt.title(ds, fontsize=25)
plt.savefig(os.path.join(img_out, ds), bbox_inches='tight')
plt.close('all')