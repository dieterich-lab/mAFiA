import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align.sam import AlignmentWriter

import matplotlib.pyplot as plt

ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch1_w_splint.fasta'
query_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_splint_batch1_m6A_RTA/basecalled.fasta'
sam_file = '/home/adrian/img_out/test.sam'
recon_ref_file = '/home/adrian/img_out/recon_ref.fasta'
out_hist_path = '/home/adrian/img_out/hist_norm_global_scores.png'

references = list(SeqIO.parse(ref_file, 'fasta'))
queries = list(SeqIO.parse(query_file, 'fasta'))

min_segment_len = 10

#########################################
### local aligner #######################
local_aligner = PairwiseAligner()
local_aligner.mode = 'local'
local_aligner.match_score = 1
local_aligner.mismatch_score = -1
local_aligner.open_gap_score = -1
local_aligner.extend_gap_score = -1

### global aligner ######################
global_aligner = PairwiseAligner()
global_aligner.mode = 'global'
global_aligner.match_score = 1
global_aligner.mismatch_score = -0.5
global_aligner.open_gap_score = -2
global_aligner.extend_gap_score = -1
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

for ind, ref in enumerate(references):
    ref.id = 'block1_{}'.format(ind)
    ref.name = ''
    ref.description = ''
    references[ind] = ref
dict_recon_references = {ref.id : ref for ref in references}

all_alignments = []
# query = queries[0]
for query in tqdm(queries):
    if len(query.seq)<min_segment_len:
        continue
    all_identified_segments = []
    remaining_seq_start = [(query.seq, 0)]

    try:
        while len(remaining_seq_start)>0:
            remaining_seq_start = remaining_seq_start[:-1] + get_daughter_seq_pos(remaining_seq_start[-1], all_identified_segments)
    except:
        continue

    ### sort by start position ###
    all_identified_segments.sort(key=lambda x: x[3])
    filtered_segments = [seg for seg in all_identified_segments if (seg[3]-seg[2])>=min_segment_len]
    # filtered_segments = all_identified_segments.copy()
    if (len(filtered_segments)==0) or (len(filtered_segments)>10):
        continue

    ### reconstruct full reference ###
    segment_sequence = [seg[0] for seg in filtered_segments]
    ref_id = 'block{}_{}'.format(len(segment_sequence), ''.join([str(s) for s in segment_sequence]))
    ref_seq = Seq('').join([references[ind].seq for ind in segment_sequence])
    ref_recon = SeqRecord(
        seq=ref_seq,
        id=ref_id,
        description=''
    )
    if ref_id not in dict_recon_references.keys():
        dict_recon_references[ref_id] = ref_recon

    ### global alignment ###
    g_align = global_aligner.align(ref_recon, query)[0]

    ### normalized global alignment score, NOT conventional mapq!!! ###
    g_align.mapq = int(g_align.score / len(query.seq) * 100)
    # print('\n=====================================================================')
    # print('Segment sequence: {}'.format('-'.join([str(i) for i in segment_sequence])))
    # print(format_alignment(*global_align, full_sequences=True))

    all_alignments.append(g_align)

### plot histogram of mapping scores ###
thresh_mapq = 30
all_mapping_scores = np.array([a.mapq for a in all_alignments])
num_pass = np.sum(all_mapping_scores>=thresh_mapq)
pass_rate = int(num_pass / len(all_mapping_scores) * 100)
plt.figure(figsize=(5, 5))
plt.hist(all_mapping_scores, range=[0, 100], bins=100)
plt.xlabel('Norm. Mapping scores', fontsize=12)
plt.ylabel('Counts', fontsize=12)
plt.axvline(x=thresh_mapq, c='g')
plt.title('Pass rate at Q$\geq${}\n{}/{} = {}%'.format(thresh_mapq, num_pass, len(all_mapping_scores), pass_rate), fontsize=15)
plt.savefig(out_hist_path, bbox_inches='tight')
plt.close('all')

### output sam file ###
with open(sam_file, 'w') as h:
    alignment_writer = AlignmentWriter(h)
    # alignment_writer.write_header(all_alignments)
    alignment_writer.write_multiple_alignments(all_alignments)

### output recon ref ###
sorted_recon_references = [dict_recon_references[k] for k in sorted(dict_recon_references.keys())]
with open(recon_ref_file, "w") as handle:
  SeqIO.write(sorted_recon_references, handle, "fasta")

# TODO for SAM output:
# CS tag
# calcs /home/adrian/img_out/test.sam -r /home/adrian/img_out/recon_ref.fasta > /home/adrian/img_out/test_cs.sam
# SAM header
