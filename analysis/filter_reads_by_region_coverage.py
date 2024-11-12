import os
import pysam
from tqdm import tqdm

ds = 'TAC'
cond = 'TAC'

chrom = 'MT'
chromStart = 3706
chromStop = 3774
strand = '+'
thresh_overlap = 69

flag_required = 0 if strand == '+' else 16

out_count = 0

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
in_bam_file = os.path.join(res_dir, f'{ds}/{cond}_merged/mt-Ti/relabelled.mAFiA.reads.bam')
out_bam_file = os.path.join(res_dir, f'{ds}/{cond}_merged/mt-Ti/thresh_overlap{thresh_overlap}.relabelled.mAFiA.reads.bam')
with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
    with pysam.AlignmentFile(out_bam_file, 'wb', template=in_bam) as out_bam:
        for this_read in tqdm(in_bam.fetch(chrom, chromStart, chromStop)):
            if this_read.flag == flag_required:
                aligned_nts = [
                    this_pair for this_pair in this_read.get_aligned_pairs(matches_only=True)
                    if (this_pair[1] >= chromStart) and (this_pair[1] <= chromStop)
                ]
                if len(aligned_nts) >= thresh_overlap:
                    out_bam.write(this_read)
                    out_count += 1
pysam.index(out_bam_file)
print(f'Written out {out_count} reads')
