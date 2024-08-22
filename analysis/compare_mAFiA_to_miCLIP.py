import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

bed6_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

miCLIP_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/miCLIP'
mAFiA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'

img_out = '/home/adrian/img_out/mAFiA_miCLIP_comparison'
os.makedirs(img_out, exist_ok=True)

cond = 'SHAM'
cond = 'TAC'

miCLIP_file = os.path.join(miCLIP_dir, f'nochr_predicted_m6A_{cond}.bed')
mAFiA_file = os.path.join(mAFiA_dir, f'{cond}_merged', 'chrALL.mAFiA.sites.bed')

df_miCLIP = pd.read_csv(miCLIP_file, sep='\t', names=bed6_fields)
df_miCLIP['name'] = 'm6A'
df_mAFiA = pd.read_csv(mAFiA_file, sep='\t', dtype={'chrom': str})
df_mAFiA.drop(columns=['score'], inplace=True)

df_merged = pd.merge(df_miCLIP, df_mAFiA, on=[
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'strand'
])

print(f'Enough coverage for {len(df_merged)} / {len(df_miCLIP)} miCLIP sites')

thresh_confidence = 80.0
df_merged_thresh = df_merged[df_merged['confidence'] >= thresh_confidence]

plt.figure(figsize=(5, 5))
plt.scatter(df_merged_thresh['score'], df_merged_thresh['modRatio'], s=1)

bin_min = 0.5
bin_max = 1.0
bin_width = 0.1
bin_edges = np.arange(bin_min, bin_max+0.01, bin_width)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

binned_modRatios = []
for bin_ind in range(len(bin_edges)-1):
    bin_start = bin_edges[bin_ind]
    bin_end = bin_edges[bin_ind+1] + 0.000001
    binned_modRatios.append(
        df_merged_thresh[(df_merged_thresh['score'] >= bin_start) * (df_merged_thresh['score'] < bin_end)]['modRatio'].values)

plt.figure(figsize=(5, 5))
plt.boxplot(binned_modRatios, positions=bin_centers, widths=0.05, whis=0.5, showfliers=False)
plt.xlim([bin_min, bin_max])
plt.xticks(bin_edges, np.round(bin_edges, 1))
plt.xlabel('miCLIP score', fontsize=12)
plt.ylabel('mAFiA stoichiometry', fontsize=12)
plt.title(cond, fontsize=15)
plt.savefig(os.path.join(img_out, f'boxplot_{cond}.png'), bbox_inches='tight')

### compare miCLIP diff. sites to mAFiA KS sites ###
