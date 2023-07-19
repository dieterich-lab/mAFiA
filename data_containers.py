import os
import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from glob import glob
from ont_fast5_api.fast5_interface import get_fast5_file
import random
random.seed(0)
from random import sample
from tqdm import tqdm

class Nucleotide:
    def __init__(self, read_id='', read_pos=-1, ref_pos=-1, pred_5mer='NNNNN', ref_5mer='NNNNN', feature=[], mod_prob=-1):
        self.read_id = str(read_id)
        self.read_pos = int(read_pos)
        self.ref_pos = int(ref_pos)
        self.pred_5mer = str(pred_5mer)
        self.ref_5mer = str(ref_5mer)
        self.feature = np.array(feature)
        self.mod_prob = float(mod_prob)

class Aligned_Read:
    def __init__(self, read_id='', read_pos=-1, ref_pos=-1, query_5mer='NNNNN', pred_5mer='NNNNN', norm_signal=[], flag=-1):
        self.read_id = str(read_id)
        self.read_pos = int(read_pos)
        self.ref_pos = int(ref_pos)
        self.query_5mer = str(query_5mer)
        self.pred_5mer = str(pred_5mer)
        self.norm_signal = np.array(norm_signal)
        self.flag = int(flag)

    def create_nucleotide(self, in_pred_5mer, in_feature):
        return Nucleotide(
            read_id = self.read_id,
            read_pos = self.read_pos,
            ref_pos = self.ref_pos,
            pred_5mer = in_pred_5mer,
            feature = in_feature,
        )

# class mRNA_Site:
#     def __init__(self, row, ref):
#         self.ind = row['index']
#         self.chr = row['Chr'].lstrip('chr')
#         self.start = row['Sites'] - 1  # 0-based
#         self.strand = row['Strand']
#         self.glori_ratio = row['Ratio']
#
#         ref_5mer = ref[self.chr][self.start-2:self.start+3]
#         if self.strand == '-':
#             self.ref_motif = str(Seq(ref_5mer).reverse_complement())
#         else:
#             self.ref_motif = ref_5mer
#
#     def print(self):
#         print('site{}, chr{}, start{}, strand{}'.format(self.ind, self.chr, self.start, self.strand))
#         print('Reference motif {}'.format(self.ref_motif))
#         print('GLORI ratio {}'.format(self.glori_ratio))

class mRNA_Site:
    def __init__(self, row, ref):
        self.chr = row['chrom']
        self.start = row['chromStart']  # 0-based
        self.strand = row['strand']
        self.ind = f'{self.chr}.{self.start}'

        ref_5mer = ref[self.chr][self.start-2:self.start+3]
        if self.strand == '-':
            self.ref_motif = str(Seq(ref_5mer).reverse_complement())
        else:
            self.ref_motif = ref_5mer

    def print(self):
        print('chr{}, start{}, strand{}'.format(self.chr, self.start, self.strand))
        print('Reference motif {}'.format(self.ref_motif))

class Data_Container:
    def __init__(self, name, bam_path, fast5_dir):
        print('Loading data {}'.format(name))
        self.name = name
        self.bam = pysam.AlignmentFile(bam_path, 'rb')
        self.nucleotides = {}

    def _index_fast5_files(self, fast5_dir, index_bam_queries_only=False):
        print('Indexing fast5 files from {}'.format(fast5_dir))
        f5_paths = glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)
        self.indexed_read_ids = {}
        if index_bam_queries_only:
            bam_query_names = [alignment.query_name for alignment in self.bam.fetch()]
        else:
            bam_query_names = []
        for f5_filepath in tqdm(f5_paths):
            try:
                f5 = get_fast5_file(f5_filepath, mode="r")
            except:
                print('Error reading {}!'.format(f5_filepath))
            else:
                read_ids = f5.get_read_ids()
                if len(bam_query_names) > 0:
                    for read_id in read_ids:
                        if read_id in bam_query_names:
                            self.indexed_read_ids[read_id] = f5_filepath
                else:
                    for read_id in read_ids:
                        self.indexed_read_ids[read_id] = f5_filepath

        print('{} reads indexed'.format(len(self.indexed_read_ids)))

    def _med_mad(self, x, factor=1.4826):
        med = np.median(x)
        mad = np.median(np.absolute(x - med)) * factor
        return med, mad

    def _get_norm_signal_from_read_id(self, id, index_paths):
        filepath = index_paths[id]
        f5 = get_fast5_file(filepath, mode="r")
        read = f5.get_read(id)
        signal = read.get_raw_data(scale=True)
        signal_start = 0
        signal_end = len(signal)
        med, mad = self._med_mad(signal[signal_start:signal_end])
        return (signal[signal_start:signal_end] - med) / mad

    def build_dict_read_ref(self):
        print('Building dictionary of reads to mapped references')
        self.dict_read_ref = {}
        for read in self.bam.fetch():
            self.dict_read_ref[read.query_name] = read.reference_name

    def flush_nts_to_dataframe(self):
        dfs = []
        for this_motif in self.nucleotides.keys():
            for nt in self.nucleotides[this_motif]:
                dfs.append(
                    pd.DataFrame(
                        [(nt.read_id, self.dict_read_ref[nt.read_id], nt.read_pos, nt.ref_pos, nt.ref_5mer, nt.pred_5mer, round(nt.mod_prob, 3))],
                        columns=['read_id', 'contig', 'read_pos', 'ref_pos', 'ref_motif', 'pred_motif', 'mod_prob']
                    )
                )
        self.nucleotides.clear()
        return pd.concat(dfs).reset_index(drop=True)

class Oligo_Data_Container(Data_Container):
    def __init__(self, name, bam_path, fast5_dir):
        super().__init__(name, bam_path, fast5_dir)
        self._index_fast5_files(fast5_dir, index_bam_queries_only=True)

    def collect_features_from_reads(self, extractor, max_num_reads):
        print('Now extracting features from {}'.format(self.name))
        if max_num_reads > 0:
            sample_read_ids = {id: self.indexed_read_ids[id] for id in
                              sample(list(self.indexed_read_ids.keys()), min(len(self.indexed_read_ids.keys()), max_num_reads))}
        else:
            sample_read_ids = self.indexed_read_ids

        read_bases_features = {}
        for query_name in tqdm(sample_read_ids.keys()):
            this_read_signal = self._get_norm_signal_from_read_id(query_name, sample_read_ids)
            this_read_features, this_read_bases = extractor.get_features_from_signal(this_read_signal)
            read_bases_features[query_name] = (this_read_bases, this_read_features)
        self.read_bases_features = read_bases_features

    def collect_motif_nucleotides(self, reference_motif, reference_generator, enforce_ref_5mer=False):
        print('Collecting nucleotides for motif {}'.format(reference_motif))

        motif_relevant_ligation_ref_ids_and_positions = reference_generator.get_motif_relevant_ligation_ref_ids_and_positions(reference_motif)
        relevant_contigs = [k for k in motif_relevant_ligation_ref_ids_and_positions.keys() if k in self.bam.references]

        this_motif_nts = []
        for contig in relevant_contigs:
            for pos in motif_relevant_ligation_ref_ids_and_positions[contig]:
                this_tPos_nts = self.collect_nucleotides_aligned_to_target_pos(contig, pos, reference_motif, enforce_ref_5mer)
                if len(this_tPos_nts) > 0:
                    this_motif_nts.extend(this_tPos_nts)

        self.nucleotides[reference_motif] = this_motif_nts
        print('{} NTs collected'.format(len(self.nucleotides[reference_motif])))

    def collect_nucleotides_aligned_to_target_pos(self, contig, target_pos, ref_motif=None, enforce_ref_5mer=False):
        all_nts = []
        for pileupcolumn in self.bam.pileup(contig, target_pos, target_pos + 1, truncate=True):
            if pileupcolumn.reference_pos == target_pos:
                valid_counts = 0
                for ind, pileupread in enumerate(pileupcolumn.pileups):
                    flag = pileupread.alignment.flag
                    if flag != 0:
                        continue
                    query_name = pileupread.alignment.query_name
                    query_position = pileupread.query_position
                    if query_position is None:
                        continue
                    query_motif = pileupread.alignment.query_sequence[(query_position - 2):(query_position + 3)]
                    if enforce_ref_5mer and (query_motif != ref_motif):
                        continue
                    if query_name in self.read_bases_features.keys():
                        this_read_bases, this_read_feature = self.read_bases_features[query_name]
                        this_site_motif = this_read_bases[(query_position - 2):(query_position + 3)]
                        if this_site_motif != query_motif:
                            print('!!! Error: Site motif {} =/= query {}!!!'.format(this_site_motif, query_motif))
                            print('Flag {}'.format(flag))
                            continue
                        this_site_feature = this_read_feature[query_position]
                        all_nts.append(Nucleotide(
                            read_id=query_name,
                            read_pos=query_position,
                            ref_pos=target_pos,
                            pred_5mer=this_site_motif,
                            ref_5mer=ref_motif,
                            feature=this_site_feature)
                        )
                        valid_counts += 1
        return all_nts

class mRNA_Data_Container(Data_Container):
    def __init__(self, name, bam_path, fast5_dir):
        super().__init__(name, bam_path, fast5_dir)
        self._index_fast5_files(fast5_dir, index_bam_queries_only=False)

    def collect_nucleotides_aligned_to_mRNA_site(self, extractor, site, thresh_coverage=1, max_num_reads=1000, enforce_ref_5mer=False):
        all_aligned_reads = []
        for pileupcolumn in self.bam.pileup(site.chr, site.start, site.start + 1, truncate=True):
            if pileupcolumn.reference_pos == site.start:
                coverage = pileupcolumn.get_num_aligned()
                if coverage > thresh_coverage:
                    valid_counts = 0
                    for pileupread in pileupcolumn.pileups:
                        flag = pileupread.alignment.flag
                        if not (
                                ((site.strand == '+') and (flag == 0))
                                or ((site.strand == '-') and (flag == 16))
                        ):
                            continue
                        query_name = pileupread.alignment.query_name
                        query_position = pileupread.query_position
                        if query_position is None:
                            continue
                        if flag == 16:
                            query_position = pileupread.alignment.query_length - query_position - 1
                        query_sequence = pileupread.alignment.get_forward_sequence()
                        query_5mer = query_sequence[(query_position - 2):(query_position + 3)]
                        if enforce_ref_5mer and (query_5mer != site.ref_motif):
                            continue
                        if (query_name in self.indexed_read_ids.keys()):
                            valid_counts += 1
                            this_read_signal = self._get_norm_signal_from_read_id(query_name, self.indexed_read_ids)
                            all_aligned_reads.append(Aligned_Read(
                                read_id=query_name,
                                read_pos=query_position,
                                query_5mer=query_5mer,
                                ref_pos=pileupcolumn.reference_pos,
                                norm_signal=this_read_signal,
                                flag=flag
                            ))
                        if (max_num_reads > 0) and (len(all_aligned_reads) >= max_num_reads):
                            break
        site_nts = extractor.get_nucleotides_from_multiple_reads(all_aligned_reads)
        for nt in site_nts:
            nt.ref_5mer = site.ref_motif
        if len(site_nts)>0:
            self.nucleotides[site.ind] = site_nts
