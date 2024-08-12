import os
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

ds = 'HEK293/100WT'
ds_name = ds.replace('/', '_')
sam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/{ds}/genome_mapped.sam'

img_out = '/home/adrian/img_out/read_lengths_mapped_vs_unmapped'
os.makedirs(img_out, exist_ok=True)

read_flag_length = {}
with pysam.AlignmentFile(sam_file, 'r') as sam:
    for this_read in tqdm(sam.fetch()):
        read_flag_length[this_read.query_name] = (this_read.flag, this_read.query_length)

unmapped_length = [this_val[1] for this_val in read_flag_length.values() if this_val[0] == 4]
mapped_length = [this_val[1] for this_val in read_flag_length.values() if this_val[0] != 4]

plt.figure(figsize=(5, 4))
plt.hist(mapped_length, range=[0, 3000], bins=50, label='mapped', alpha=0.5)
plt.hist(unmapped_length, range=[0, 3000], bins=50, label='unmapped', alpha=0.5)
plt.legend(loc='upper right')
plt.xlabel('Read length', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.title(ds_name, fontsize=15)
plt.savefig(os.path.join(img_out, f'{ds_name}.png'))
