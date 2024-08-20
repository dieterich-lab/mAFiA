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

# thresh_conf = 50.0
# thresh_modRatio = 50.0

# thresh_conf = 0.0
# thresh_modRatio = 0.0

# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']

# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

for (thresh_conf, thresh_modRatio) in [(0.0, 0.0), (50.0, 50.0)]:
    for this_cond in conditions:
        in_bed_file = os.path.join(res_dir, ds, this_cond, 'chrALL.mAFiA.sites.bed')

        df_in = pd.read_csv(in_bed_file, sep='\t', dtype={'chrom': str})
        df_in_thresh = df_in[
            (df_in['modRatio'] >= thresh_modRatio)
            * (df_in['confidence'] >= thresh_conf)
            ]

        for this_mod in ['m6A', 'psi']:
            out_bed_file = os.path.join(res_dir, 'metaPlotR', f'{ds}_{this_cond}_{this_mod}_modRatio{thresh_modRatio}.bed')

            df_out = df_in_thresh[df_in_thresh['name'] == this_mod]
            df_out['chrom'] = ['chr' + this_chr for this_chr in df_out['chrom']]
            df_out['score'] = df_out['modRatio']

            df_out[bed6_fields].to_csv(out_bed_file, sep='\t', header=False, index=False)
