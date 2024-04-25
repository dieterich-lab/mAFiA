import os
import pandas as pd

annot_dir = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102'
in_gtf_file = os.path.join(annot_dir, 'GRCh38.102.gtf')

gtf_fields = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
]

gtf = pd.read_csv(in_gtf_file, sep='\t', names=gtf_fields, dtype={'seqname': str})

for chr in [str(i) for i in range(1, 23)] + ['X']:
    chr_gtf = gtf[gtf['seqname']==chr]
    chr_gtf.to_csv(os.path.join(annot_dir, f'chr{chr}.GRCh38.102.gtf'), header=False, index=False)
