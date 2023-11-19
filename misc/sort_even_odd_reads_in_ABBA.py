import os
import pysam
import numpy as np
import pandas as pd

workspace = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/WUE_ABBA/mAFiA'
bamfile = os.path.join(workspace, 'mAFiA.reads.bam')
bam = pysam.AlignmentFile(bamfile, 'rb')

block_size = 26
thresh_diff = 192

even_reads = []
odd_reads = []
readID_readPos_refPos_modProb = []
for read in bam.fetch():
    readID = read.query_name
    readPos_modProb = list(read.modified_bases_forward.values())[0]
    read_ref_position = {qpos : rpos for qpos, rpos in read.get_aligned_pairs()}
    refPos_modProb = [(read_ref_position[readPos], modProb) for (readPos, modProb) in readPos_modProb if read_ref_position[readPos] is not None]

    even_modProbs = [tup[1] for tup in refPos_modProb if tup[0]%block_size==6]
    odd_modProbs = [tup[1] for tup in refPos_modProb if tup[0]%block_size==19]

    if len(even_modProbs)>=4 and len(odd_modProbs)>=4:
        diff = np.mean(even_modProbs) - np.mean(odd_modProbs)
        if diff > thresh_diff:
            even_reads.append(read)
        elif diff < -thresh_diff:
            odd_reads.append(read)

    for read_tup, ref_tup in zip(readPos_modProb, refPos_modProb):
        readID_readPos_refPos_modProb.append((readID, read_tup[0], ref_tup[0], ref_tup[1]/255.0))


with pysam.AlignmentFile(os.path.join(workspace, 'even.bam'), 'wb', header=bam.header) as outf:
    for read in even_reads:
        outf.write(read)

with pysam.AlignmentFile(os.path.join(workspace, 'odd.bam'), 'wb', header=bam.header) as outf:
    for read in odd_reads:
        outf.write(read)

df_parsed = pd.DataFrame(readID_readPos_refPos_modProb, columns=['read_id', 'read_pos', 'ref_pos', 'mod_prob'])
df_parsed.to_csv(os.path.join(workspace, 'pos_modProb.tsv'), sep='\t', index=False)