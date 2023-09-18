import os, re
import numpy as np
import pandas as pd
from calcs import trim, annotate, call_cstag
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, Alignment
from Bio.Align.sam import AlignmentWriter

class OligoReferenceGenerator:
    def __init__(self, oligo_ref_file=None, ligation_ref_file=None, annotation_file=None, mod_type=None, oligo_fmt='-M([0-9]+)S([0-9]+)'):
        self.oligo_fmt = oligo_fmt
        if ligation_ref_file:
            self.ligation_ref = {}
            for ligation_ref in list(SeqIO.parse(ligation_ref_file, 'fasta')):
                ligation_ref.name = ''
                ligation_ref.description = ''
                self.ligation_ref[ligation_ref.id] = ligation_ref
            self.oligo_ref = [
                self.ligation_ref[k]
                for k, v in self.ligation_ref.items()
                if len(re.findall(self.oligo_fmt, v.id))==1
            ]
        else:
            self.oligo_ref = list(SeqIO.parse(oligo_ref_file, 'fasta'))
            self.ligation_ref = {}
            for idx, ref in enumerate(self.oligo_ref.copy()):
                ref.id = self._generate_ligation_ref_id([idx])
                ref.name = ''
                ref.description = ''
                self.ligation_ref[ref.id] = ref
        self.oligo_lens = [len(ref.seq) for ref in self.oligo_ref]
        self.min_segment_len = min([len(ref.seq) for ref in self.oligo_ref]) // 2

        if annotation_file:
            self.annotation = pd.read_csv(annotation_file, sep='\t')
            if mod_type=='mod':
                self.annotation = self.annotation[self.annotation['modRatio']==100]
            elif mod_type=='unm':
                self.annotation = self.annotation[self.annotation['modRatio']==0]

    def _generate_ligation_ref_id(self, seg_seq):
        oligo_ids = [self.oligo_ref[i].id for i in seg_seq]
        oligo_origins = [id.split('-')[0] for id in oligo_ids]
        if len(set(oligo_origins))==1:
            oligo_origin = oligo_origins[0]
        else:
            raise Exception('Error - multiple origins for oligos')
        oligo_suffixes = [id.split('-')[1] for id in oligo_ids]
        return '{}-{}'.format(oligo_origin, '-'.join(oligo_suffixes))

    def get_ligation_reference(self, segments):
        segment_sequence = [seg[0] for seg in segments]
        ligation_ref_id = self._generate_ligation_ref_id(segment_sequence)
        ligation_ref_seq = Seq('').join([self.oligo_ref[ind].seq for ind in segment_sequence])
        ligation_ref_record = SeqRecord(
            seq=ligation_ref_seq,
            id=ligation_ref_id,
            description=''
        )
        if ligation_ref_id not in self.ligation_ref.keys():
            self.ligation_ref[ligation_ref_id] = ligation_ref_record

        return ligation_ref_record

    def collect_motif_oligos(self):
        motif_oligo_name_size_pos = {}

        if self.annotation is None:
            for this_oligo_ref in self.oligo_ref:
                motif_seq_indices = re.findall(self.oligo_fmt, this_oligo_ref.id)
                if len(motif_seq_indices) == 1:
                    this_oligo_name = this_oligo_ref.id
                    # this_oligo_index = motif_seq_indices[0][0]
                    this_oligo_size = len(this_oligo_ref.seq)
                    this_oligo_location = this_oligo_size // 2
                    this_oligo_motif = str(this_oligo_ref.seq[this_oligo_location-2 : this_oligo_location+3])
                    if this_oligo_motif not in motif_oligo_name_size_pos.keys():
                        motif_oligo_name_size_pos[this_oligo_motif] = {}
                    motif_oligo_name_size_pos[this_oligo_motif][this_oligo_name] = (this_oligo_size, this_oligo_location)
        else:
            for _, this_row in self.annotation.iterrows():
                this_oligo_name, this_oligo_location, this_modRatio, this_oligo_motif = this_row[['chrom', 'chromStart', 'modRatio', 'ref5mer']]
                this_oligo_size = len(self.ligation_ref[this_oligo_name].seq)
                if this_oligo_motif not in motif_oligo_name_size_pos.keys():
                    motif_oligo_name_size_pos[this_oligo_motif] = {}
                motif_oligo_name_size_pos[this_oligo_motif][this_oligo_name] = (this_oligo_size, this_oligo_location)


        self.motif_oligos = motif_oligo_name_size_pos
        self.flat_oligo_dims = {kk: vv for k, v in self.motif_oligos.items() for kk, vv in v.items()}

    def get_motif_relevant_ligation_ref_ids_and_positions(self, in_motif):
        relevant_oligos = self.motif_oligos[in_motif]
        relevant_ligation_ref_ids_positions = {}
        for ref_id in self.ligation_ref.keys():
            for this_oligo_id in relevant_oligos.keys():
                oligo_origin, oligo_suffix = this_oligo_id.split('-')
                lig_ref_decomposed = ref_id.split('-')
                lig_ref_origin = lig_ref_decomposed[0]
                lig_ref_suffices = lig_ref_decomposed[1:]
                if (lig_ref_origin==oligo_origin) and (oligo_suffix in lig_ref_suffices):
                    site_positions = []
                    for lig_ind, this_lig_ref_suffix in enumerate(lig_ref_suffices):
                        if this_lig_ref_suffix==oligo_suffix:
                            if lig_ind>0:
                                previous_suffices = lig_ref_suffices[:lig_ind] if lig_ind>0 else 0
                                inter_block_shifts = sum([self.flat_oligo_dims['{}-{}'.format(lig_ref_origin, block)][0] for block in previous_suffices])
                            else:
                                inter_block_shifts = 0
                            intra_block_shift = self.flat_oligo_dims['{}-{}'.format(lig_ref_origin, this_lig_ref_suffix)][1]
                            site_positions.append(inter_block_shifts+intra_block_shift)
                    relevant_ligation_ref_ids_positions[ref_id] = site_positions

        return relevant_ligation_ref_ids_positions

class QueryContainer:
    def __init__(self, query_file):
        format = os.path.basename(query_file).split('.')[-1]
        self.records = list(SeqIO.parse(query_file, format=format))
        if format == 'fastq':
            for query in self.records:
                query.letter_annotations = {}

    def get_records(self):
        return self.records

    def __len__(self):
        return len(self.records)

class LocalAligner(PairwiseAligner):
    def __init__(self):
        super().__init__()
        self.mode = 'local'
        self.match_score = 1
        self.mismatch_score = -1
        self.open_gap_score = -1
        self.extend_gap_score = -1

    def get_best_local_segment(self, in_seq, ref_generator):
        ref_alignments = [
            self.align(in_seq, ref.seq)[0]
            for ref in ref_generator.oligo_ref
        ]
        ref_scores = [a.score for a in ref_alignments]
        chosen_ref_ind = np.argmax(ref_scores)
        chosen_score = ref_scores[chosen_ref_ind]
        chosen_alignment = ref_alignments[chosen_ref_ind]
        target_start = chosen_alignment.coordinates[0][0]
        target_end = chosen_alignment.coordinates[0][-1]
        query_start = chosen_alignment.coordinates[1][0]
        query_end = chosen_alignment.coordinates[1][-1]

        return chosen_ref_ind, chosen_score, target_start, target_end, chosen_alignment

class Global_Aligner(PairwiseAligner):
    def __init__(self):
        super().__init__()
        self.mode = 'global'
        self.match_score = 1
        self.mismatch_score = -1
        self.open_gap_score = -5
        self.extend_gap_score = -1

    def get_recon_align_by_global_alignment(self, query, ref_recon):
        ### global alignment ###
        recon_align = self.align(ref_recon, query)[0]

        ### normalized global alignment score, NOT conventional mapq!!! ###
        recon_align.mapq = int(recon_align.score / len(query.seq) * 100)
        # print('\n=====================================================================')
        # print('Segment sequence: {}'.format('-'.join([str(i) for i in segment_sequence])))
        # print(format_alignment(*global_align, full_sequences=True))

        return recon_align

class Splitter:
    def __init__(self, in_aligner, min_segment_len=10, thresh_mapq=50, homopolymer=False):
        self.aligner = in_aligner
        self.min_segment_len = int(min_segment_len)
        self.thresh_mapq = int(thresh_mapq)
        self.homopolymer = bool(homopolymer)
        if self.homopolymer:
            print('Homopolymer only', flush=True)

    def _get_daughter_seq_pos(self, mother_seq_pos, identified_segments, ref_generator, min_fragment_len=8):
        mother_seq, mother_pos = mother_seq_pos
        daughter_seq_pos = []
        best_ref_ind, best_score, best_start, best_end, best_alignment = self.aligner.get_best_local_segment(mother_seq,
                                                                                                        ref_generator)
        identified_segments.append(
            (best_ref_ind, best_score, best_start + mother_pos, best_end + mother_pos, best_alignment))
        if best_start >= min_fragment_len:
            daughter_seq_pos.append((mother_seq[:best_start], mother_pos))
        if len(mother_seq) - best_end >= min_fragment_len:
            daughter_seq_pos.append((mother_seq[best_end:], mother_pos + best_end))
        return daughter_seq_pos

    def split_query_into_segments(self, in_query, ref_generator):
        all_identified_segments = []
        remaining_seq_start = [(in_query.seq, 0)]

        try:
            while len(remaining_seq_start) > 0:
                remaining_seq_start = remaining_seq_start[:-1] + self._get_daughter_seq_pos(remaining_seq_start[-1],
                                                                                      all_identified_segments,
                                                                                      ref_generator)
                # if args.debug:
                #     print('Identified segments:', all_identified_segments)
                #     print('Remaining segments', remaining_seq_start)
        except:
            return []

        ### sort segments by start position ###
        all_identified_segments.sort(key=lambda x: x[3])
        out_filtered_segments = [seg for seg in all_identified_segments
                             if (((seg[3] - seg[2]) >= self.min_segment_len) and (
                        (seg[1] / (seg[3] - seg[2]) * 100) >= self.thresh_mapq))]
        # out_filtered_segments = all_identified_segments.copy()
        if (len(out_filtered_segments) == 0) or (len(out_filtered_segments) > 10):
            return []
        # if args.debug:
        #     print('All identified segments:', all_identified_segments)
        #     for seg in out_filtered_segments:
        #         print('Target {} - {}'.format(seg[2], seg[3]))
        #         print(seg[-1])

        ### check homopolymer ###
        if self.homopolymer and (len(np.unique([seg[0] for seg in out_filtered_segments])) > 1):
            return []

        return out_filtered_segments

class Chainer:
    def get_recon_align_by_chain(self, in_segments, in_target, in_query):
        full_query_seq = in_query.seq

        ### fill gaps ###
        padded_segments = []
        curr_pos = 0
        # curr_pos = in_segments[0][2]
        while len(in_segments) > 0:
            seg = in_segments[0]
            begin, end = seg[2:4]
            if curr_pos < begin:
                padded_segments.append(
                    (None, None, curr_pos, begin, None)
                )
                curr_pos = begin
            else:
                padded_segments.append(seg)
                curr_pos = end
                in_segments = in_segments[1:]
        if curr_pos < len(full_query_seq):
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
                recon_target_lines.append('-' * (q_end - q_start))

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

        recon_align.mapq = int(np.mean([seg[1] / (seg[3] - seg[2]) for seg in padded_segments if seg[1]]) * 100)
        recon_align.target = in_target
        recon_align.query = in_query

        return recon_align

class Writer:
    def __init__(self, out_sam_file, out_ligation_ref_file, write_md, write_cs):
        self.out_sam_file = out_sam_file
        self.out_ligation_ref_file = out_ligation_ref_file
        self.write_md = write_md
        self.write_cs = write_cs

    def get_correct_sam_line(self, in_alignment, sam_writer, write_md=False, write_cs=True):
        sam_fields = sam_writer.format_alignment(in_alignment, md=write_md).rstrip('\n').split('\t')

        cigar_string = sam_fields[5]
        match_in_del_sclip = re.findall(r"([0-9]+)M|([0-9]+)I|([0-9]+)D|([0-9]+)S", cigar_string)
        if len(match_in_del_sclip[0][1]) > 0:
            match_in_del_sclip[0] = ('', '', '', match_in_del_sclip[0][1])
            if len(match_in_del_sclip[1][2]) > 0:
                match_in_del_sclip = [match_in_del_sclip[0]] + match_in_del_sclip[2:]
        if len(match_in_del_sclip[-1][1]) > 0:
            match_in_del_sclip[-1] = ('', '', '', match_in_del_sclip[-1][1])
            if len(match_in_del_sclip[-2][2]) > 0:
                match_in_del_sclip = match_in_del_sclip[:-2] + [match_in_del_sclip[-1]]
        if len(match_in_del_sclip[0][2]) > 0:
            match_in_del_sclip = match_in_del_sclip[1:]
        if len(match_in_del_sclip[-1][2]) > 0:
            match_in_del_sclip = match_in_del_sclip[:-1]

        correct_cigar_string = ''.join(
            [n + suffix for cs in match_in_del_sclip for n, suffix in zip(cs, ['M', 'I', 'D', 'S']) if len(n) > 0])
        sam_fields[5] = correct_cigar_string

        str_match_stick = in_alignment._format_unicode().split('\n')[1]
        correct_pos = in_alignment.indices[0, str_match_stick.find('|')] + 1  # 1-based
        sam_fields[3] = str(correct_pos)

        if write_cs:
            len_clips = trim.get_softclip_lengths(correct_cigar_string)
            query_trimmed = trim.softclips(sam_fields[9], len_clips)
            ref_trimmed = trim.unmapped_region(str(in_alignment.target.seq), correct_pos - 1, correct_cigar_string)
            ref_anno = annotate.insertion(ref_trimmed, correct_cigar_string)
            query_anno = annotate.deletion(query_trimmed, correct_cigar_string)
            cs_tag = call_cstag.short_form(call_cstag.long_form(ref_anno, query_anno))
            sam_fields += [cs_tag]

        ### reassemble ###
        out_sam_line = '\t'.join(sam_fields) + '\n'

        return out_sam_line

    def _sort_ref_ids(self, ref_ids, ref_generator):
        return [
            ref_generator.ligation_ref[k] for k in sorted(ref_generator.ligation_ref.keys(),
                                                          key=lambda x: (len(x.split('-')), x.split('-')[1])
                                                          ) if k in ref_ids
        ]

    def write_sam(self, alignments, ref_ids, ref_generator):
        sorted_recon_references = self._sort_ref_ids(ref_ids, ref_generator)

        print('Writing to samfile...', flush=True)
        with open(self.out_sam_file, 'w') as out_sam:
            ### write header SQ lines ###
            for ref in sorted_recon_references:
                out_sam.write('@SQ\tSN:{}\tLN:{}\n'.format(ref.id, len(ref.seq)))
            out_sam.write('@PG\tID:spomlette\tPN:spomlette\tVN:0.01\tCL:blah\n')

            ### write read alignments ###
            alignment_writer = AlignmentWriter(out_sam)

            for this_alignment in alignments:
                corrected_sam_line = self.get_correct_sam_line(this_alignment, alignment_writer,
                                                                 write_md=self.write_md, write_cs=self.write_cs)
                # corrected_sam_line = alignment_writer.format_alignment(this_alignment)
                out_sam.write(corrected_sam_line)

        ### output recon ref ###
        with open(self.out_ligation_ref_file, "w") as handle:
            SeqIO.write(sorted_recon_references, handle, "fasta")
