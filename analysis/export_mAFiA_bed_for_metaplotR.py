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

# polyA = ['below_90', 'above_90']

# ds = 'HFpEF'
# conditions = ['ctrl', 'HFpEF']

# ds = 'Diet'
# conditions = ['WT_CD', 'WT_WD']

# ds = 'TAC'
# conditions = ['SHAM', 'TAC']

ds = 'CM'
conditions = ['WT', 'M3KO']

in_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
out_dir = os.path.join(in_dir, 'metaPlotR')

for this_cond in conditions:
    in_bed_file = os.path.join(in_dir, ds, f'{this_cond}_merged', 'chrALL.mAFiA.sites.bed')
    for (thresh_conf, thresh_modRatio) in [(0.0, 0.0), (50.0, 50.0)]:
        df_in = pd.read_csv(in_bed_file, sep='\t', dtype={'chrom': str})
        df_in_thresh = df_in[
            (df_in['modRatio'] >= thresh_modRatio)
            * (df_in['confidence'] >= thresh_conf)
            ]

        for this_mod in ['m6A', 'psi']:
            out_bed_file = os.path.join(in_dir, 'metaPlotR',
                                        f'{ds}_{this_cond}_{this_mod}_modRatio{thresh_modRatio}.bed')
            df_out = df_in_thresh[df_in_thresh['name'] == this_mod].copy()
            df_out['chrom'] = ['chr' + this_chr for this_chr in df_out['chrom']]
            df_out['score'] = df_out['modRatio']

            df_out[bed6_fields].to_csv(out_bed_file, sep='\t', header=False, index=False)
