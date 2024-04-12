import os
import pysam

workspace = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT/chr1'
in_bam_file = os.path.join(workspace, 'chr1_100082344_100082361.mAFiA.reads.bam')

sort_pos = 100082357
mod_code = 17802
out_bam_high_file = os.path.join(workspace, f'chr1_100082344_100082361.highAt{sort_pos}.mAFiA.reads.bam')
out_bam_low_file = os.path.join(workspace, f'chr1_100082344_100082361.lowAt{sort_pos}.mAFiA.reads.bam')

reads = {}
with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
    for this_read in in_bam.fetch():
        this_read.set_tag('MM', this_read.get_tag('MM').replace('N+21891', 'N+a'))
        reads[this_read.query_name] = this_read

id_scores = []
for this_id, this_read in reads.items():
    # print(this_read.query_name, this_read.modified_bases)
    dict_read_to_ref_pos = {tup[0]: tup[1] for tup in this_read.get_aligned_pairs(matches_only=True)}
    for this_mod_base in this_read.modified_bases[('N', 0, mod_code)]:
        pos = dict_read_to_ref_pos[this_mod_base[0]]
        if pos==sort_pos:
            score = this_mod_base[1] / 255.0
            id_scores.append((this_id, score))
            break
# id_scores.sort(key=lambda x: x[1], reverse=True)

with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
    with pysam.AlignmentFile(out_bam_high_file, 'wb', template=in_bam) as out_bam_high:
        with pysam.AlignmentFile(out_bam_low_file, 'wb', template=in_bam) as out_bam_low:
            for id, score in id_scores:
                if score>=0.5:
                    out_bam_high.write(reads[id])
                else:
                    out_bam_low.write(reads[id])

pysam.index(out_bam_high_file)
pysam.index(out_bam_low_file)