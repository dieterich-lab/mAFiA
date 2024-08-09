import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os


benchmark = 'PRAISE'
benchmark_file = '/home/adrian/Data/PRAISE/PRAISE_HEK293T_span1-3.bed'
wt_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/chrALL.mAFiA.sites.bed'
oe_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293_TRUB1_OE/merged/chrALL.mAFiA.sites.bed'
kd_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/NanoSPA/HEK_siTRUB1_input_merged/chrALL.mAFiA.sites.bed'

img_out = '/home/adrian/img_out/psi-co-mAFiA'

df_benchmark = pd.read_csv(benchmark_file, sep='\t', dtype={'chrom': str})
df_wt = pd.read_csv(wt_file, sep='\t', dtype={'chrom': str}).drop(columns=['score']).rename(
    columns={'coverage': 'coverage_WT', 'modRatio': 'modRatio_WT', 'confidence': 'confidence_WT'})
df_oe = pd.read_csv(oe_file, sep='\t', dtype={'chrom': str}).drop(columns=['score']).rename(
    columns={'coverage': 'coverage_OE', 'modRatio': 'modRatio_OE', 'confidence': 'confidence_OE'})
df_kd = pd.read_csv(kd_file, sep='\t', dtype={'chrom': str}).drop(columns=['score']).rename(
    columns={'coverage': 'coverage_KD', 'modRatio': 'modRatio_KD', 'confidence': 'confidence_KD'})

TRUB1_motifs = ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT']
df_benchmark_TRUB1 = df_benchmark[df_benchmark['ref5mer'].isin(TRUB1_motifs)]

df_merged = pd.merge(df_benchmark_TRUB1, df_wt, how='outer', on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'])
df_merged = pd.merge(df_merged, df_oe, how='outer', on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'])
df_merged = pd.merge(df_merged, df_kd, how='outer', on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ref5mer'])
df_merged = df_merged[df_merged['ref5mer'].isin(TRUB1_motifs)]

num_sites = {
    'benchmark': (~df_merged['score'].isna()).sum(),
    'benchmark-WT': ((~df_merged['score'].isna()) * (~df_merged['modRatio_WT'].isna())).sum(),
    'benchmark-OE': ((~df_merged['score'].isna()) * (~df_merged['modRatio_OE'].isna())).sum(),
    'benchmark-KD': ((~df_merged['score'].isna()) * (~df_merged['modRatio_KD'].isna())).sum(),
    'WT': (~df_merged['modRatio_WT'].isna()).sum(),
    'WT-OE': ((~df_merged['modRatio_WT'].isna()) * (~df_merged['modRatio_OE'].isna())).sum(),
    'WT-KD': ((~df_merged['modRatio_WT'].isna()) * (~df_merged['modRatio_KD'].isna())).sum()
}

condition_names = {
    'WT': 'HEK293-WT',
    'OE': 'TRUB1-OE',
    'KD': 'TRUB1-KD'
}

xylim = [-5, 105]
plt.figure(figsize=(15, 10))
for subplot_ind, cond in enumerate(['WT', 'OE', 'KD']):
    plt.subplot(2, 3, subplot_ind+1)
    plt.plot(df_merged['score'], df_merged[f'modRatio_{cond}'], '.')
    plt.plot([0, 100], [0, 100], 'r--')
    plt.xlim(xylim)
    plt.ylim(xylim)
    plt.xlabel(rf'$S_{{{benchmark}}}$', fontsize=12)
    plt.ylabel(rf'$S_{{{condition_names[cond]}}}$', fontsize=12)
    plt.title(f"{num_sites[f'benchmark-{cond}']} / {num_sites['benchmark']}")
for subplot_ind, cond in enumerate(['OE', 'KD']):
    plt.subplot(2, 3, 3+subplot_ind+2)
    plt.plot(df_merged[f'modRatio_WT'], df_merged[f'modRatio_{cond}'], '.')
    plt.plot([0, 100], [0, 100], 'r--')
    plt.xlim(xylim)
    plt.ylim(xylim)
    plt.xlabel(rf"$S_{{{'WT'}}}$", fontsize=12)
    plt.ylabel(rf'$S_{{{condition_names[cond]}}}$', fontsize=12)
    plt.title(f"{num_sites[f'WT-{cond}']} / {num_sites['WT']}")
plt.suptitle('GUUCN\nAll sites with overage$\geq$20', fontsize=15)
plt.savefig(os.path.join(img_out, 'comparison_PRAISE_WT_OE_KD.png'), bbox_inches='tight')