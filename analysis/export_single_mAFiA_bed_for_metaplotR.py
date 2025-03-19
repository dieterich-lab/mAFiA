import os
import pandas as pd

bed6_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'

ds = 'WT'
in_bed_file = os.path.join(base_dir, 'HEK293/WT_P2/chrALL.mAFiA.sites.bed')
out_dir = os.path.join(base_dir, 'HEK293/WT_P2/metaPlotR')

os.makedirs(out_dir, exist_ok=True)

for (thresh_conf, thresh_modRatio) in [(0.0, 0.0), (50.0, 50.0)]:
    df_in = pd.read_csv(in_bed_file, sep='\t', dtype={'chrom': str})
    df_in_thresh = df_in[
        (df_in['modRatio'] >= thresh_modRatio)
        * (df_in['confidence'] >= thresh_conf)
        ]

    for this_mod in ['m6A', 'psi']:
        out_bed_file = os.path.join(out_dir,
                                    f'{ds}_{this_mod}_modRatio{thresh_modRatio}.bed')
        df_out = df_in_thresh[df_in_thresh['name'] == this_mod].copy()
        df_out['chrom'] = ['chr' + this_chr for this_chr in df_out['chrom']]
        df_out['score'] = df_out['modRatio']

        df_out[bed6_fields].to_csv(out_bed_file, sep='\t', header=False, index=False)
