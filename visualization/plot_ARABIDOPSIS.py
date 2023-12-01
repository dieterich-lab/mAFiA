import pandas as pd

mAFiA_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Arabidopsis_thaliana/col0/chr1/mAFiA.sites.bed'
miclip_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/parker_miclip_sites.tds'

df_mAFiA = pd.read_csv(mAFiA_file, sep='\t', dtype={'chrom': str})
df_miclip = pd.read_csv(miclip_file, sep='\t',
                        usecols=range(3),
                        names=['chrom', 'chromStart', 'chromEnd'],
                        dtype={'chrom': str})

df_merged = pd.merge(df_mAFiA, df_miclip, on=['chrom', 'chromStart', 'chromEnd'])