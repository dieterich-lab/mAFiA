import os
import numpy as np
import pysam
from glob import glob
from utils import index_fast5_files
import random
random.seed(0)
from random import sample
from tqdm import tqdm
from utils import get_norm_signal_from_read_id

class nucleotide:
    def __init__(self, read_id='', read_pos=-1, ref_pos=-1, pred_5mer='NNNNN', ref_5mer='NNNNN', feature=[], mod_prob=-1):
        self.read_id = str(read_id)
        self.read_pos = int(read_pos)
        self.ref_pos = int(ref_pos)
        self.pred_5mer = str(pred_5mer)
        self.ref_5mer = str(ref_5mer)
        self.feature = np.array(feature)
        self.mod_prob = float(mod_prob)

class aligned_read:
    def __init__(self, read_id='', read_pos=-1, ref_pos=-1, query_5mer='NNNNN', pred_5mer='NNNNN', norm_signal=[], flag=-1):
        self.read_id = str(read_id)
        self.read_pos = int(read_pos)
        self.ref_pos = int(ref_pos)
        self.query_5mer = str(query_5mer)
        self.pred_5mer = str(pred_5mer)
        self.norm_signal = np.array(norm_signal)
        self.flag = int(flag)

    def create_nucleotide(self, in_pred_5mer, in_feature):
        return nucleotide(
            read_id = self.read_id,
            read_pos = self.read_pos,
            ref_pos = self.ref_pos,
            pred_5mer = in_pred_5mer,
            feature = in_feature,
        )

class data_container:
    def __init__(self, bam_path, fast5_dir):
        self.bam = pysam.AlignmentFile(bam_path, 'rb')
        self.f5_paths = glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)
        print('Indexing fast5 files from {}'.format(fast5_dir))
        self.indexed_read_ids = index_fast5_files(self.f5_paths, self.bam)
        print('{} reads indexed'.format(len(self.indexed_read_ids)))

class oligo_data_container(data_container):
    def __init__(self, bam_path, fast5_dir):
        super().__init__(bam_path, fast5_dir)
        self.motif_nts = {}

    def collect_features_from_reads(self, extractor, max_num_reads):
        if max_num_reads > 0:
            sample_read_ids = {id: self.indexed_read_ids[id] for id in
                              sample(list(self.indexed_read_ids.keys()), min(len(self.indexed_read_ids.keys()), max_num_reads))}
        else:
            sample_read_ids = self.indexed_read_ids

        read_bases_features = {}
        for query_name in tqdm(sample_read_ids.keys()):
            this_read_signal = get_norm_signal_from_read_id(query_name, sample_read_ids)
            this_read_features, this_read_bases = extractor.get_features_from_signal(this_read_signal)
            read_bases_features[query_name] = (this_read_bases, this_read_features)
        self.read_bases_features = read_bases_features

    def collect_motif_nucleotides(self, motif_ind, ref, block_size, block_center, enforce_ref_5mer=False):
        relevant_contigs = [k for k in ref.keys() if ((motif_ind in k.split('_')[1]) and (k in self.bam.references))]
        this_motif_nts = []
        for contig in relevant_contigs:
            block_str = contig.split('_')[1]
            site_positions = np.where(np.array(list(block_str)) == motif_ind)[0] * block_size + block_center
            for pos in site_positions:
                reference_motif = ref[contig][(pos - 2):(pos + 3)]
                this_tPos_nts = self.collect_nucleotides_aligned_to_target_pos(contig, pos, reference_motif, enforce_ref_5mer)
                if len(this_tPos_nts) > 0:
                    this_motif_nts.extend(this_tPos_nts)

        self.motif_nts[motif_ind] = this_motif_nts
        print('{} NTs collected'.format(len(self.motif_nts[motif_ind])))

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
                        all_nts.append(nucleotide(
                            read_id=query_name,
                            read_pos=query_position,
                            ref_pos=target_pos,
                            pred_5mer=this_site_motif,
                            ref_5mer=ref_motif,
                            feature=this_site_feature)
                        )
                        valid_counts += 1
        return all_nts

class mRNA_data_container(data_container):
    def __init__(self, bam_path, fast5_dir):
        super().__init__(bam_path, fast5_dir)