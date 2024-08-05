import argparse
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np

# infile = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/chrALL.mAFiA.sites.bed'
# outfile = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/chrALL.bedMethyl'

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile')
parser.add_argument('-o', '--outfile')
args = parser.parse_args()
print(args)

mAFiA_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
    'coverage',
    'modRatio',
    'confidence'
]

bedMethyl_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'thickStart',
    'thickEnd',
    'itemRgb',
    'coverage',
    'frequency'
]

dict_mod_names = {
    'm6A': 'm6A',
    'psi': 'Y'
}

in_bed = pd.read_csv(args.infile, sep='\t', dtype={'chrom': str})

out_bed = in_bed.copy()
out_bed['name'] = [dict_mod_names[this_name] for this_name in in_bed['name']]
log_pval = -np.log10((100.0 - in_bed['confidence']) / 100.0)
log_pval[np.isinf(log_pval)] = 1000
out_bed['score'] = np.int64(np.round(log_pval))
out_bed['thickStart'] = out_bed['chromStart']
out_bed['thickEnd'] = out_bed['chromEnd']
out_bed['itemRgb'] = '0,0,0'
out_bed['frequency'] = np.int64(np.round(in_bed['modRatio']))

out_bed = out_bed[bedMethyl_fields]
out_bed.to_csv(args.outfile, sep='\t', index=False)