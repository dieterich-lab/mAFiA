# requires:
# pip install biopython==1.81
# pip install calcs

import numpy as np
import argparse
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, Alignment
from Bio.Align.sam import AlignmentWriter
import os, re
import matplotlib.pyplot as plt
from calcs import trim, annotate, call_cstag

parser = argparse.ArgumentParser(description='Map oligo basecalls through the Spanish omelette method')
parser.add_argument("--ref_file", type=str)
parser.add_argument("--query_file", type=str)
parser.add_argument("--recon_ref_file", type=str)
parser.add_argument("--sam_file", type=str)
parser.add_argument("--thresh_mapq", type=int, default=50)
parser.add_argument("--homopolymer", default=False, action="store_true")
parser.add_argument("--debug", default=False, action="store_true")
parser.add_argument("--write_md", default=False, action="store_true")
parser.add_argument("--write_cs", default=False, action="store_true")
args = parser.parse_args()

THRESH_MAPQ = int(args.thresh_mapq)

# ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch1_w_splint.fasta'
# query_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_splint_batch1_m6A_RTA/basecalled.fasta'
# recon_ref_file = '/home/adrian/img_out/spomlette_recon_ref_q{}.fasta'.format(THRESH_MAPQ)
# sam_file = '/home/adrian/img_out/spomlette_q{}.sam'.format(THRESH_MAPQ)
# out_hist_path = '/home/adrian/img_out/hist_norm_global_scores_q{}.png'.format(THRESH_MAPQ)
out_hist_path = os.path.join(os.path.dirname(args.sam_file), 'hist_spomlette_mapping_scores_q{}.png'.format(THRESH_MAPQ))

references = list(SeqIO.parse(args.ref_file, 'fasta'))
queries = list(SeqIO.parse(args.query_file, 'fasta'))

BLOCK_SIZE = np.min([len(ref.seq) for ref in references])
MIN_SEGMENT_LEN = BLOCK_SIZE // 2

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
global_aligner.mismatch_score = -1
global_aligner.open_gap_score = -5
global_aligner.extend_gap_score = -1
#########################################

def get_local_segment(in_seq):
    ref_alignments = [
        local_aligner.align(in_seq, ref.seq)[0]
        for ref in references
    ]
    ref_scores = [a.score for a in ref_alignments]
    chosen_ref_ind = np.argmax(ref_scores)
    chosen_score = ref_scores[chosen_ref_ind]
    chosen_alignment = ref_alignments[chosen_ref_ind]
    target_start = chosen_alignment.coordinates[0][0]
    target_end = chosen_alignment.coordinates[0][-1]
    query_start = chosen_alignment.coordinates[1][0]
    query_end = chosen_alignment.coordinates[1][-1]

    ### debug ###
    # if args.debug:
    #     chosen_seq = chosen_alignment.query.seq
    #     chosen_motif = chosen_seq[len(chosen_seq)//2-2:len(chosen_seq)//2+3]
    #     print('Match scores: {}'.format(', '.join([str(score) for score in ref_scores])))
    #     print('Chosen reference {} - {}'.format(chosen_ref_ind, chosen_motif))
    #     print(chosen_alignment)

    return chosen_ref_ind, chosen_score, target_start, target_end, chosen_alignment

def get_daughter_seq_pos(mother_seq_pos, identified_segments, min_fragment_len=8):
    mother_seq, mother_pos = mother_seq_pos
    daughter_seq_pos = []
    best_ref_ind, best_score, best_start, best_end, best_alignment = get_local_segment(mother_seq)
    identified_segments.append((best_ref_ind, best_score, best_start+mother_pos, best_end+mother_pos, best_alignment))
    if best_start >= min_fragment_len:
        daughter_seq_pos.append((mother_seq[:best_start], mother_pos))
    if len(mother_seq)-best_end >= min_fragment_len:
        daughter_seq_pos.append((mother_seq[best_end:], mother_pos+best_end))
    return daughter_seq_pos

def get_recon_align_by_global_alignment(in_filtered_segments):
    ### global alignment ###
    recon_align = global_aligner.align(ref_recon, query)[0]

    ### normalized global alignment score, NOT conventional mapq!!! ###
    recon_align.mapq = int(recon_align.score / len(query.seq) * 100)
    # print('\n=====================================================================')
    # print('Segment sequence: {}'.format('-'.join([str(i) for i in segment_sequence])))
    # print(format_alignment(*global_align, full_sequences=True))

    return recon_align

def get_recon_align_by_chain(in_segments, in_target, in_query):
    full_query_seq = in_query.seq

    ### fill gaps ###
    padded_segments = []
    curr_pos = 0
    # curr_pos = in_segments[0][2]
    while len(in_segments)>0:
        seg = in_segments[0]
        begin, end = seg[2:4]
        if curr_pos<begin:
            padded_segments.append(
                (None, None, curr_pos, begin, None)
            )
            curr_pos = begin
        else:
            padded_segments.append(seg)
            curr_pos = end
            in_segments = in_segments[1:]
    if curr_pos<len(full_query_seq):
        padded_segments.append(
            (None, None, curr_pos, len(full_query_seq), None)
        )

    ### chain together segments ###
    recon_query_lines = []
    recon_target_lines = []
    for seg in padded_segments:
        ref_ind, local_score, q_start, q_end = seg[:4]
        seg_alignment = seg[-1]

        if seg_alignment is not None:
            t_start = seg_alignment.coordinates[1][0]
            t_end = seg_alignment.coordinates[1][-1]
            recon_query_lines.append(
                ''.join(['-' for i in range(t_start)]) \
                + str(seg_alignment[0]) \
                + ''.join(['-' for i in range(t_end, len(seg_alignment.query))])
            )

            recon_target_lines.append(
                str(seg_alignment.query[:t_start]) \
                + str(seg_alignment[1]) \
                + str(seg_alignment.query[t_end:])
            )

        else:
            recon_query_lines.append(str(full_query_seq[q_start:q_end]))
            recon_target_lines.append('-'*(q_end-q_start))

    # if args.debug:
    #     for i in range(len(padded_segments)):
    #         print('\n')
    #         print('\t'*4+' '*3, recon_target_lines[i])
    #         print('\t'*4+' '*3, recon_query_lines[i])
    #         print(padded_segments[i][-1])

    recon_query_aligned = ''.join(recon_query_lines)
    recon_target_aligned = ''.join(recon_target_lines)

    ### create alignment object ###
    recon_lines = [recon_target_aligned, recon_query_aligned]

    ### full coords ###
    recon_coords = Alignment.infer_coordinates(recon_lines)
    recon_seqs = [l.replace('-', '') for l in recon_lines]
    recon_align = Alignment(recon_seqs, recon_coords)

    ### reduced coords to generate correct soft-clipping ###
    # reduced_coords = recon_coords[:, (recon_coords[0] != 0) * (recon_coords[0] < np.max(recon_coords[0]))]
    # recon_align = Alignment(recon_seqs, reduced_coords)

    recon_align.mapq = int(np.mean([seg[1] / (seg[3]-seg[2]) for seg in padded_segments if seg[1]]) * 100)
    recon_align.target = in_target
    recon_align.query = in_query

    return recon_align

def get_m6A_line(alignment, ref_id, block_size):
    line_aligned = alignment._format_unicode().split('\n')[0]
    q_start = alignment.coordinates[0, 0]
    cumsum_non_gap = q_start + np.cumsum(np.int32([x!='-' for x in list(line_aligned)])) - 1

    num_blocks = int(ref_id.split('_')[0].lstrip('block'))
    m6A_pos = block_size // 2 + block_size * np.arange(num_blocks)
    m6A_line = [' '] * len(line_aligned)
    for pos in m6A_pos:
        if pos in cumsum_non_gap.flatten():
            shifted_pos = np.where(cumsum_non_gap==pos)[0][0]
            for sp in range(shifted_pos-2, shifted_pos+3):
                if (sp>=0) and (sp<len(m6A_line)):
                    m6A_line[sp] = '*'
    return ''.join(m6A_line)

def get_correct_sam_line(in_alignment, sam_writer, write_md=False, write_cs=True):
    sam_fields = sam_writer.format_alignment(in_alignment, md=write_md).rstrip('\n').split('\t')

    ### cigar string: change I on both ends to S ###
    cigar_string = sam_fields[5]
    match_in_del_sclip = re.findall(r"([0-9]+)M|([0-9]+)I|([0-9]+)D|([0-9]+)S", cigar_string)
    # beginning I #
    if len(match_in_del_sclip[0][1]) > 0:
        match_in_del_sclip[0] = ('', '', '', match_in_del_sclip[0][1])
        ### remove 'D' after starting 'S' ###
        if len(match_in_del_sclip[1][2]) > 0:
            match_in_del_sclip = [match_in_del_sclip[0]] + match_in_del_sclip[2:]
    # ending I #
    if len(match_in_del_sclip[-1][1]) > 0:
        match_in_del_sclip[-1] = ('', '', '', match_in_del_sclip[-1][1])
        ### remove 'D' before ending 'S' ###
        if len(match_in_del_sclip[-2][2]) > 0:
            match_in_del_sclip = match_in_del_sclip[:-2] + [match_in_del_sclip[-1]]
    # beginning D #
    if len(match_in_del_sclip[0][2]) > 0:
        match_in_del_sclip = match_in_del_sclip[1:]
    # ending D #
    if len(match_in_del_sclip[-1][2]) > 0:
        match_in_del_sclip = match_in_del_sclip[:-1]

    correct_cigar_string = ''.join([n+suffix for cs in match_in_del_sclip for n,suffix in zip(cs, ['M', 'I', 'D', 'S']) if len(n)>0])
    sam_fields[5] = correct_cigar_string

    ### start pos ###
    str_match_stick = in_alignment._format_unicode().split('\n')[1]
    correct_pos = in_alignment.indices[0, str_match_stick.find('|')] + 1   # 1-based
    sam_fields[3] = str(correct_pos)

    if args.debug:
        correct_MIDS = re.findall(r"([0-9]+)M|([0-9]+)I|([0-9]+)D|([0-9]+)S", correct_cigar_string)
        ### check SM...MS ###
        if (len(correct_MIDS[0][3])>0) or (len(correct_MIDS[-1][3])>0)\
                or (len(correct_MIDS[0][2])>0) or (len(correct_MIDS[-1][2])>0):
            print(in_alignment)
            print('Start pos:', correct_pos)
            print('Before:', cigar_string)
            print('After:', correct_cigar_string)

            tlen = 0
            qlen = 0
            for x in correct_MIDS:
                m, i, d, s = [int(xx) if len(xx) else 0 for xx in x]
                tlen += m + d
                qlen += m + i + s
            print('Target len:', len(in_alignment.target)==tlen)
            print('Query len:', len(in_alignment.query)==qlen)

            print('\n')

    ### generate CS flag ###
    if write_cs:
        len_clips = trim.get_softclip_lengths(correct_cigar_string)
        query_trimmed = trim.softclips(sam_fields[9], len_clips)
        ref_trimmed = trim.unmapped_region(str(in_alignment.target.seq), correct_pos-1, correct_cigar_string)
        ref_anno = annotate.insertion(ref_trimmed, correct_cigar_string)
        query_anno = annotate.deletion(query_trimmed, correct_cigar_string)
        cs_tag = call_cstag.short_form(call_cstag.long_form(ref_anno, query_anno))
        sam_fields += [cs_tag]

    ### reassemble ###
    out_sam_line = '\t'.join(sam_fields) + '\n'

    return out_sam_line

dict_recon_references = {}
dict_motifs = {}
for idx, ref in enumerate(references.copy()):
    motif_ind = int(re.findall(r"M([0-9]+)", ref.id)[0])
    dict_motifs[idx] = motif_ind
    ref.id = 'block1_{}'.format(motif_ind)
    ref.name = ''
    ref.description = ''
    dict_recon_references[ref.id] = ref
# dict_recon_references = {ref.id : ref for ref in references}

all_alignments = []
print('Now mapping reads...', flush=True)
for query in tqdm(queries):
    if len(query.seq)<MIN_SEGMENT_LEN:
        continue
    all_identified_segments = []
    remaining_seq_start = [(query.seq, 0)]

    try:
        while len(remaining_seq_start)>0:
            remaining_seq_start = remaining_seq_start[:-1] + get_daughter_seq_pos(remaining_seq_start[-1], all_identified_segments)
            # if args.debug:
            #     print('Identified segments:', all_identified_segments)
            #     print('Remaining segments', remaining_seq_start)
    except:
        continue

    ### sort segments by start position ###
    all_identified_segments.sort(key=lambda x: x[3])
    filtered_segments = [seg for seg in all_identified_segments
                         if (((seg[3]-seg[2])>=MIN_SEGMENT_LEN) and ((seg[1]/(seg[3]-seg[2])*100)>=THRESH_MAPQ))]
    # filtered_segments = all_identified_segments.copy()
    if (len(filtered_segments)==0) or (len(filtered_segments)>10):
        continue
    # if args.debug:
    #     print('All identified segments:', all_identified_segments)
    #     for seg in filtered_segments:
    #         print('Target {} - {}'.format(seg[2], seg[3]))
    #         print(seg[-1])

    ### check homopolymer ###
    if (args.homopolymer) and (len(np.unique([seg[0] for seg in filtered_segments]))>1):
        continue

    ### reconstruct full reference ###
    segment_sequence = [seg[0] for seg in filtered_segments]
    ref_id = 'block{}_{}'.format(len(segment_sequence), ''.join([str(dict_motifs[s]) for s in segment_sequence]))
    ref_seq = Seq('').join([references[ind].seq for ind in segment_sequence])
    ref_recon = SeqRecord(
        seq=ref_seq,
        id=ref_id,
        description=''
    )
    if ref_id not in dict_recon_references.keys():
        dict_recon_references[ref_id] = ref_recon

    ### concatenate local alignments ###
    # full_alignment = get_recon_align_by_global_alignment(filtered_segments)
    full_alignment = get_recon_align_by_chain(filtered_segments, ref_recon, query)

    if args.debug:
        m6A_line = get_m6A_line(full_alignment, ref_id, block_size=BLOCK_SIZE)
        print('\n')
        print(query.id)
        print('Inferred reference:', ref_id)
        print(m6A_line)
        print(full_alignment._format_unicode())
        print('Score {}'.format(full_alignment.mapq))

    all_alignments.append(full_alignment)

### plot histogram of mapping scores ###
all_mapping_scores = np.array([a.mapq for a in all_alignments])
num_pass = np.sum(all_mapping_scores>=THRESH_MAPQ)
pass_rate = int(num_pass / len(queries) * 100)
plt.figure(figsize=(5, 5))
plt.hist(all_mapping_scores, range=[0, 100], bins=100)
plt.xlabel('Chain Mapping scores', fontsize=12)
plt.ylabel('Counts', fontsize=12)
plt.axvline(x=THRESH_MAPQ, c='g')
plt.title('Pass rate at Q$\geq${}\n{}/{} = {}%'.format(THRESH_MAPQ, num_pass, len(queries), pass_rate), fontsize=15)
plt.savefig(out_hist_path, bbox_inches='tight')
plt.close('all')

### filter alignments by score ###
filtered_alignments = [a for a in all_alignments if a.mapq>=THRESH_MAPQ]
filtered_ref_ids = [a.target.id for a in filtered_alignments]

### output sam file ###
sorted_recon_references = [
    dict_recon_references[k] for k in sorted(dict_recon_references.keys(),
                                             key=lambda x: (int(x.split('_')[0].lstrip('block')), x.split('_')[1]))
    if k in filtered_ref_ids
]

print('Writing to samfile...', flush=True)
with open(args.sam_file, 'w') as out_sam:
    ### write header SQ lines ###
    for ref in sorted_recon_references:
        out_sam.write('@SQ\tSN:{}\tLN:{}\n'.format(ref.id, len(ref.seq)))
    out_sam.write('@PG\tID:spomlette\tPN:spomlette\tVN:0.01\tCL:blah\n')

    ### write read alignments ###
    alignment_writer = AlignmentWriter(out_sam)
    # alignment_writer.write_header(all_alignments)
    # alignment_writer.write_multiple_alignments(filtered_alignments)

    for this_alignment in filtered_alignments:
        corrected_sam_line = get_correct_sam_line(this_alignment, alignment_writer, write_md=args.write_md, write_cs=args.write_cs)
        # corrected_sam_line = alignment_writer.format_alignment(this_alignment)
        out_sam.write(corrected_sam_line)

### output recon ref ###
with open(args.recon_ref_file, "w") as handle:
  SeqIO.write(sorted_recon_references, handle, "fasta")
