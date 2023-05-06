from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from Bio.Align.sam import AlignmentWriter
import numpy as np
from tqdm import tqdm

from Bio.Align import PairwiseAligner

ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch1_w_splint.fasta'
query_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_splint_batch1_m6A_RTA/basecalled.fasta'

references = list(SeqIO.parse(ref_file, 'fasta'))
queries = list(SeqIO.parse(query_file, 'fasta'))

min_segment_len = 10

#########################################
### local aligner ###
local_aligner = PairwiseAligner()
local_aligner.mode = 'local'
local_aligner.match_score = 1
local_aligner.mismatch_score = -1
local_aligner.open_gap_score = -1
local_aligner.extend_gap_score = -1

### global aligner ######################
global_aligner = PairwiseAligner()
global_aligner.mode = 'global'
global_aligner.match_score = 2
global_aligner.mismatch_score = -1
global_aligner.open_gap_score = -5
global_aligner.extend_gap_score = -2
#########################################

def get_most_likely_segment(in_seq):
    ref_alignments = [
        local_aligner.align(in_seq, ref)[0]
        for ref in references
    ]
    ref_scores = [a.score for a in ref_alignments]
    chosen_ref_ind = np.argmax(ref_scores)
    chosen_score = ref_scores[chosen_ref_ind]
    chosen_alignment = ref_alignments[chosen_ref_ind]
    q_start = chosen_alignment.coordinates[0][0]
    q_end = chosen_alignment.coordinates[0][-1]

    return chosen_ref_ind, chosen_score, q_start, q_end

def get_daughter_seq_pos(mother_seq_pos, identified_segments, min_fragment_len=8):
    mother_seq, mother_pos = mother_seq_pos
    daughter_seq_pos = []
    best_ref_ind, best_score, best_start, best_end = get_most_likely_segment(mother_seq)
    identified_segments.append((best_ref_ind, best_score, best_start+mother_pos, best_end+mother_pos))
    if best_start >= min_fragment_len:
        daughter_seq_pos.append((mother_seq[:best_start], mother_pos))
    if len(mother_seq)-best_end >= min_fragment_len:
        daughter_seq_pos.append((mother_seq[best_end:], mother_pos+best_end))
    return daughter_seq_pos

# query = queries[0]
with open('/home/adrian/img_out/dummy.sam', 'w') as h:
    alignment_writer = AlignmentWriter(h, md=True)
    for query in tqdm(queries[:10]):
        if len(query.seq)<min_segment_len:
            continue
        all_identified_segments = []
        remaining_seq_start = [(query.seq, 0)]

        while len(remaining_seq_start)>0:
            remaining_seq_start = remaining_seq_start[:-1] + get_daughter_seq_pos(remaining_seq_start[-1], all_identified_segments)

        ### sort by start position ###
        all_identified_segments.sort(key=lambda x: x[3])
        filtered_segments = [seg for seg in all_identified_segments if (seg[3]-seg[2])>=min_segment_len]
        # filtered_segments = all_identified_segments.copy()
        if len(filtered_segments)==0:
            continue

        ### reconstruct full reference ###
        segment_sequence = [seg[0] for seg in filtered_segments]
        ref_recon = Seq('').join([references[ind].seq for ind in segment_sequence])

        ### global alignment ###
        # g_align1 = pairwise2.align.globalms(ref_recon, query.seq, 2, -1, -5, -2, one_alignment_only=True)[0]
        g_align = global_aligner.align(ref_recon, query.seq)[0]

        # print('\n=====================================================================')
        # print('Segment sequence: {}'.format('-'.join([str(i) for i in segment_sequence])))
        # print(format_alignment(*global_align, full_sequences=True))

        alignment_writer.write_single_alignment([g_align])