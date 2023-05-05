from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np

ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch1_w_splint.fasta'
query_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_splint_batch1_m6A_RTA/basecalled.fasta'

references = list(SeqIO.parse(ref_file, 'fasta'))
queries = list(SeqIO.parse(query_file, 'fasta'))

def get_most_likely_segment(in_seq, pos_shift=0):
    ref_alignments = [
        pairwise2.align.localms(in_seq, ref, 2, -1, -2, -2, one_alignment_only=True)[0]
        for ref in references
    ]
    ref_scores = [a.score for a in ref_alignments]
    chosen_ref_ind = np.argmax(ref_scores)
    chosen_score = ref_scores[chosen_ref_ind]
    chosen_alignment = ref_alignments[chosen_ref_ind]
    formatted = format_alignment(*chosen_alignment)
    print(formatted)
    format_start, format_str = [x for x in formatted.split(' ') if len(x)>0][0:2]
    q_start = int(format_start) - 1
    q_len = len(format_str.rstrip('\n').replace('-', ''))
    q_end = q_start + q_len

    return chosen_ref_ind, chosen_score, (q_start+pos_shift), (q_end+pos_shift)

def get_daughter_seq_pos(mother_seq_pos, identified_segments):
    mother_seq, mother_pos = mother_seq_pos
    daughter_seq_pos = []
    best_ref_ind, best_score, best_start, best_end = get_most_likely_segment(mother_seq)
    identified_segments.append((best_ref_ind, best_score, best_start+mother_pos, best_end+mother_pos))
    if best_start>=10:
        daughter_seq_pos.append((mother_seq[:best_start], mother_pos))
    if (len(mother_seq)-best_end)>=10:
        daughter_seq_pos.append((mother_seq[best_end:], mother_pos+best_end))
    return daughter_seq_pos

query = queries[0]
full_seq = query.seq

all_identified_segments = []
remaining_segments = [(full_seq, 0)]

while len(remaining_segments)>0:
    new_segments = get_daughter_seq_pos(remaining_segments[-1], all_identified_segments)
    remaining_segments = remaining_segments[:-1] + new_segments

### sort by start position ###
all_identified_segments.sort(key=lambda x: x[3])