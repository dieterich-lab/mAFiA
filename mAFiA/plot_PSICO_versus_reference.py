import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import Counter
import os

img_out = '/home/adrian/img_out/PSICO'

# ds = '100_WT_0_IVT'
ds = 'P2_WT'
# ref = 'BID-Seq'
ref = 'PRAISE'
res_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/psU/PSICO_inference/HEK293/{ds}/chosen8_{ref}/psUco.sites.bed'

df_res = pd.read_csv(res_file, sep='\t')

# df_exact = pd.read_csv('/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/psU/site_annotations/BID-Seq_PRAISE_merged.bed', sep='\t')
# df_res_exact = pd.merge(df_res, df_exact, on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'])

# motif_counts = Counter(df_res['ref5mer']).most_common()

motifs = [
    'GGTGG',
    'GTTCA',
    'GTTCC',
    'TGTGG',
    'TGTAG',
    'AGTGG',
    'GTTCG',
    'GGTCC'
]

num_motifs = 8
num_rows = 2
num_cols = 4

min_coverage = 20

plt.figure(figsize=(16, 8))
for subplot_ind, motif in enumerate(motifs):
    # sub_df = df_res_exact[
    #     (df_res_exact['ref5mer']==motif)
    #     * (df_res_exact['coverage']>=min_coverage)
    #     ]
    sub_df = df_res[
        (df_res['ref5mer']==motif)
        * (df_res['coverage']>=min_coverage)
    ]
    plt.subplot(num_rows, num_cols, subplot_ind+1)
    plt.plot(sub_df['score'], sub_df['modRatio'], '.')
    plt.plot([-5, 105], [-5, 105], 'r--', alpha=0.5)
    plt.xlim([-5, 105])
    plt.ylim([-5, 105])
    plt.title(motif)
    if subplot_ind>=(num_cols*(num_rows-1)):
        plt.xlabel(ref, fontsize=12)
    if subplot_ind%num_cols==0:
        plt.ylabel('$\psi$-mAFiA', fontsize=12)
plt.suptitle(f'{ds}\n{ref}', fontsize=15)
plt.savefig(os.path.join(img_out, f'scatter_PSICO_vs_{ref}_{ds}_{num_motifs}motifs_{min_coverage}cov.png'), bbox_inches='tight')
plt.close('all')