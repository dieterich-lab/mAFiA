import os
import argparse
import pysam
import numpy as np
import re

# prj_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian'
# # dataset = 'A_RTA'
# dataset = 'm6A_RTA'
# infile = os.path.join(prj_dir, '{}_sorted.bam'.format(dataset))
# outfile = infile.replace('sorted.bam', 'sorted_filtered.bam')
# indel_thresh = 10

parser = argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--outfile')
parser.add_argument('--indel_thresh', default=10)
args = parser.parse_args()
indel_thresh = int(args.indel_thresh)
print('Filtering {} with max indel per read {}'.format(args.infile, indel_thresh))

### check bam file ###
bam_in = pysam.AlignmentFile(args.infile, 'rb')
bam_out = pysam.AlignmentFile(args.outfile, "wb", template=bam_in)

in_counts = 0
out_counts = 0
thresholded_reads = []
for read in bam_in.fetch():
    in_counts += 1
    if read.flag in [0, 2048]:
        cs = read.cigarstring
        res = re.split('(\d+)', cs)[1:]
        cs_pairs = list(zip(np.int32(res[0::2]), res[1::2]))
        ins_len = [pair[0] for pair in cs_pairs if pair[1]=='I']
        del_len = [pair[0] for pair in cs_pairs if pair[1]=='D']
        max_ins_len = np.max(ins_len) if (len(ins_len)>0) else 0
        max_del_len = np.max(del_len) if (len(del_len)>0) else 0

        if (max_ins_len<indel_thresh) and (max_del_len<indel_thresh):
            out_counts += 1
            bam_out.write(read)
bam_in.close()
bam_out.close()

survival_rate = out_counts / in_counts
print('Before: {} reads'.format(in_counts))
print('After: {} reads'.format(out_counts))
print('Survival rate: {:.2f}'.format(survival_rate))