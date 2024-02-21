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
import h5py
from joblib import Parallel, delayed

class Nucleotide:
    def __init__(self, read_id='', read_pos=-1, ref_pos=-1, pred_5mer='NNNNN', ref_5mer='NNNNN', feature=None, strand='.', mod_type='N', mod_prob=-1):
        self.read_id = str(read_id)
        self.read_pos = int(read_pos)
        self.ref_pos = int(ref_pos)
        self.pred_5mer = str(pred_5mer)
        self.ref_5mer = str(ref_5mer)
        if feature is None:
            self.feature = np.array([])
        else:
            self.feature = np.array(feature)
        self.strand = str(strand)
        self.mod_type = str(mod_type)
        self.mod_prob = float(mod_prob)

class Aligned_Read:
    def __init__(self, read_id='', read_pos=-1, ref_pos=-1, query_5mer='NNNNN', pred_5mer='NNNNN', norm_signal=None, flag=-1, strand='.'):
        self.read_id = str(read_id)
        self.read_pos = int(read_pos)
        self.ref_pos = int(ref_pos)
        self.query_5mer = str(query_5mer)
        self.pred_5mer = str(pred_5mer)
        if norm_signal is None:
            self.norm_signal = np.array([])
        else:
            self.norm_signal = np.array(norm_signal)
        self.flag = int(flag)
        self.strand = str(strand)

    def create_nucleotide(self, in_pred_5mer, in_feature):
        return Nucleotide(
            read_id = self.read_id,
            read_pos = self.read_pos,
            ref_pos = self.ref_pos,
            pred_5mer = in_pred_5mer,
            feature = in_feature,
            strand = self.strand
        )

class mRNASite:
    def __init__(self, row, ref):
        self.chr = row['chrom']
        self.start = row['chromStart']  # 0-based
        self.strand = row['strand']
        self.ind = f'{self.chr}.{self.start}'

        ref_5mer = ref[self.chr][self.start-2:self.start+3]
        if self.strand == '-':
            self.ref_5mer = str(Seq(ref_5mer).reverse_complement())
        else:
            self.ref_5mer = ref_5mer

    def print(self):
        print(f'chr{self.chr}, start{self.start}, strand{self.strand}', flush=True)
        print(f'Reference motif {self.ref_5mer}', flush=True)

class DataContainer:
    def __init__(self, name, bam_path):
        print(f'Loading data {name}')
        self.name = name
        self.bam = pysam.AlignmentFile(bam_path, 'rb')
        self.nucleotides = {}

    def _index_fast5_files(self, fast5_dir, index_bam_queries_only=False):
        print(f'Indexing fast5 files from {fast5_dir}')
        f5_paths = glob(os.path.join(fast5_dir, '*.fast5'), recursive=True)
        self.indexed_read_ids = {}
        if index_bam_queries_only:
            bam_query_names = [alignment.query_name for alignment in self.bam.fetch()]
        else:
            bam_query_names = []
        for f5_filepath in tqdm(f5_paths):
            try:
                with get_fast5_file(f5_filepath, mode="r") as f5:
                    read_ids = f5.get_read_ids()
                    if len(bam_query_names) > 0:
                        for read_id in read_ids:
                            if read_id in bam_query_names:
                                self.indexed_read_ids[read_id] = f5_filepath
                    else:
                        for read_id in read_ids:
                            self.indexed_read_ids[read_id] = f5_filepath
            except Exception as e:
                print(e)
                print(f'Error reading {f5_filepath}!')

        print(f'{len(self.indexed_read_ids)} reads indexed')

    def _med_mad(self, x, factor=1.4826):
        med = np.median(x)
        mad = np.median(np.absolute(x - med)) * factor
        return med, mad

    def _get_norm_signal_from_read_id(self, id, index_paths):
        filepath = index_paths[id]
        with get_fast5_file(filepath, mode="r") as f5:
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

class OligoDataContainer(DataContainer):
    def __init__(self, name, bam_path, fast5_dir):
        super().__init__(name, bam_path)
        self._index_fast5_files(fast5_dir, index_bam_queries_only=True)

    def collect_features_from_reads(self, extractor, max_num_reads):
        print(f'Now extracting features from {self.name}')
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
        print(f'Collecting nucleotides for motif {reference_motif}')

        motif_relevant_ligation_ref_ids_and_positions = reference_generator.get_motif_relevant_ligation_ref_ids_and_positions(reference_motif)
        relevant_contigs = [k for k in motif_relevant_ligation_ref_ids_and_positions.keys() if k in self.bam.references]

        this_motif_nts = []
        for contig in relevant_contigs:
            for pos in motif_relevant_ligation_ref_ids_and_positions[contig]:
                this_tPos_nts = self.collect_nucleotides_aligned_to_target_pos(contig, pos, reference_motif, enforce_ref_5mer)
                if len(this_tPos_nts) > 0:
                    this_motif_nts.extend(this_tPos_nts)

        self.nucleotides[reference_motif] = this_motif_nts
        print(f'{len(self.nucleotides[reference_motif])} NTs collected')

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
                            # print(f'!!! Error: Site motif {this_site_motif} =/= query {query_motif}!!!')
                            # print(f'Flag {flag}')
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


class MultiReadContainer(DataContainer):
    def __init__(self, name, in_bam_path, fast5_dir):
        print('Starting with fast5')
        super().__init__(name, in_bam_path)
        self._index_fast5_files(fast5_dir, index_bam_queries_only=False)


    def collect_nucleotides_on_single_read(self, read, read_features, df_sites):
        all_nts = []
        query_name = read.query_name
        flag = read.flag
        read_len = len(read.seq)
        if flag==0:
            dict_ref_to_read_pos = {tup[1]: tup[0] for tup in read.get_aligned_pairs() if (tup[0] is not None) and (tup[1] is not None)}
            matching_strand = '+'
        elif flag==16:
            dict_ref_to_read_pos = {tup[1]: (read_len-tup[0]-1) for tup in read.get_aligned_pairs() if (tup[0] is not None) and (tup[1] is not None)}
            matching_strand = '-'
        else:
            return all_nts

        sub_df_sites = df_sites[
            (df_sites['chrom'] == read.reference_name)
            * (df_sites['chromStart'] >= (read.reference_start+2))
            * (df_sites['chromEnd'] < (read.reference_end-2))
            * df_sites['chromStart'].isin(dict_ref_to_read_pos.keys())
            * df_sites['strand'] == matching_strand
            ]

        try:
            for _, row in sub_df_sites.iterrows():
                chromStart = row['chromStart']
                strand = row['strand']
                ref5mer = row['ref5mer']
                mod_type = row['name']

                query_pos = dict_ref_to_read_pos[chromStart]
                query_5mer = read.get_forward_sequence()[query_pos-2:query_pos+3]
                this_site_feature = read_features[query_pos]
                all_nts.append(Nucleotide(
                    read_id=query_name,
                    read_pos=query_pos,
                    ref_pos=chromStart,
                    strand=strand,
                    pred_5mer=query_5mer,
                    ref_5mer=ref5mer,
                    feature=this_site_feature,
                    mod_type=mod_type
                ))
        except:
            print(f'Error in {read.query_name}', flush=True)
            pass

        return all_nts

    def _get_matching_nucleotide_from_row(self, in_row, in_read, in_dict, in_read_features):
        chromStart = in_row['chromStart']
        strand = in_row['strand']
        ref5mer = in_row['ref5mer']
        mod_type = in_row['name']

        query_pos = in_dict[chromStart]
        query_5mer = in_read.get_forward_sequence()[query_pos - 2:query_pos + 3]
        this_site_feature = in_read_features[query_pos]

        return Nucleotide(
            read_id=in_read.query_name,
            read_pos=query_pos,
            ref_pos=chromStart,
            strand=strand,
            pred_5mer=query_5mer,
            ref_5mer=ref5mer,
            feature=this_site_feature,
            mod_type=mod_type
        )

    def parallel_collect_nucleotides_on_single_read(self, read, read_features, df_sites, num_jobs=16):
        all_nts = []
        flag = read.flag
        read_len = len(read.seq)
        if flag==0:
            dict_ref_to_read_pos = {tup[1]: tup[0] for tup in read.get_aligned_pairs() if ((tup[0] is not None) and (tup[1] is not None))}
            matching_strand = '+'
        elif flag==16:
            dict_ref_to_read_pos = {tup[1]: (read_len-tup[0]-1) for tup in read.get_aligned_pairs() if ((tup[0] is not None) and (tup[1] is not None))}
            matching_strand = '-'
        else:
            return all_nts

        dict_ref_to_read_pos = {k: v for k, v in dict_ref_to_read_pos.items() if v < read_len}

        sub_df_sites = df_sites[
            (df_sites['chrom'] == read.reference_name)
            * (df_sites['chromStart'] >= (read.reference_start+2))
            * (df_sites['chromEnd'] < (read.reference_end-2))
            * df_sites['chromStart'].isin(dict_ref_to_read_pos.keys())
            * df_sites['strand']==matching_strand
            ]
        num_jobs = min(len(sub_df_sites), num_jobs)
        all_nts = Parallel(n_jobs=num_jobs)(delayed(self._get_matching_nucleotide_from_row)(sub_df_sites.iloc[i], read, dict_ref_to_read_pos, read_features) for i in range(len(sub_df_sites)))

        return all_nts

    def process_reads(self, extractor, df_sites, multimod_motif_classifiers, sam_writer, write_chunk_size=100):
        if os.path.exists(sam_writer.out_sam_path):
            processed_reads = sam_writer.get_processed_reads()
        else:
            processed_reads = []
        processed_read_ids = []
        sam_writer.open()
        if len(processed_reads)>0:
            for this_read in processed_reads:
                sam_writer.fo.write(this_read)
                processed_read_ids.append(this_read.query_name)
            print(f'Skipping {len(processed_read_ids)} reads')

        reads_mod_nts = []
        for this_read in tqdm(self.bam.fetch()):
            if this_read.query_name in processed_read_ids:
                continue
            if this_read.flag not in [0, 16]:
                continue
            this_read_signal = self._get_norm_signal_from_read_id(this_read.query_name, self.indexed_read_ids)
            this_read_features, this_read_bases = extractor.get_features_from_signal(this_read_signal)
            if len(this_read_bases)!=len(this_read.seq):
                continue
            this_read_nts = self.collect_nucleotides_on_single_read(this_read, this_read_features, df_sites)
            # this_read_nts = self.parallel_collect_nucleotides_on_single_read(this_read, this_read_features, df_sites)

            mod_motif_nts = {}
            for this_nt in this_read_nts:
                if this_nt.mod_type not in mod_motif_nts.keys():
                    mod_motif_nts[this_nt.mod_type] = {}
                if this_nt.ref_5mer not in mod_motif_nts[this_nt.mod_type].keys():
                    mod_motif_nts[this_nt.mod_type][this_nt.ref_5mer] = [this_nt]
                else:
                    mod_motif_nts[this_nt.mod_type][this_nt.ref_5mer].append(this_nt)

            out_mod_nts = {}
            for this_mod in mod_motif_nts.keys():
                this_mod_nts = []
                for this_motif, nts in mod_motif_nts[this_mod].items():
                    nt_features = np.vstack([this_nt.feature for this_nt in nts])
                    mod_probs = multimod_motif_classifiers[this_mod][this_motif].binary_model.predict_proba(nt_features)[:, 1]
                    for this_nt, this_mod_prob in zip(nts, mod_probs):
                        this_nt.mod_prob = this_mod_prob
                        this_mod_nts.append(this_nt)
                out_mod_nts[this_mod] = this_mod_nts
            reads_mod_nts.append((this_read, out_mod_nts))
            if len(reads_mod_nts)>=write_chunk_size:
                sam_writer.write_reads(reads_mod_nts)
                reads_mod_nts = []
        sam_writer.write_reads(reads_mod_nts)
        sam_writer.close()


    def _get_mod_prob_nt(self, this_nt, multimod_motif_classifiers):
        this_nt.mod_prob = multimod_motif_classifiers[this_nt.mod_type][this_nt.ref_5mer].binary_model.predict_proba(this_nt.feature)[0, 1]
        return this_nt

    def process_reads_parallel(self, extractor, df_sites, multimod_motif_classifiers, sam_writer, num_jobs=16):
        for this_read in tqdm(self.bam.fetch()):
            if this_read.flag not in [0, 16]:
                continue
            this_read_signal = self._get_norm_signal_from_read_id(this_read.query_name, self.indexed_read_ids)
            this_read_features, this_read_bases = extractor.get_features_from_signal(this_read_signal)
            this_read_nts = self.collect_nucleotides_on_single_read(this_read, this_read_features, df_sites)

            mod_prob_nts = Parallel(n_jobs=num_jobs)(delayed(self._get_mod_prob_nt)(this_read_nts[i], multimod_motif_classifiers) for i in range(len(this_read_nts)))
            out_mod_nts = {}
            for this_nt in mod_prob_nts:
                if this_nt.mod_type not in out_mod_nts.keys():
                    out_mod_nts[this_nt.mod_type] = {}
                if this_nt.ref_5mer not in out_mod_nts[this_nt.mod_type].keys():
                    out_mod_nts[this_nt.mod_type][this_nt.ref_5mer] = [this_nt]
                else:
                    out_mod_nts[this_nt.mod_type][this_nt.ref_5mer].append(this_nt)

            sam_writer.write_read(this_read, out_mod_nts)


class mRNADataContainer(DataContainer):
    def __init__(self, name, bam_path, fast5_dir):
        print('Starting with fast5')
        super().__init__(name, bam_path)
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
                        if enforce_ref_5mer and (query_5mer != site.ref_5mer):
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
                                flag=flag,
                                strand=site.strand
                            ))
                        if (max_num_reads > 0) and (len(all_aligned_reads) >= max_num_reads):
                            break
        site_nts = extractor.get_nucleotides_from_multiple_reads(all_aligned_reads)
        for nt in site_nts:
            nt.ref_5mer = site.ref_5mer
        if len(site_nts)>0:
            self.nucleotides[site.ind] = site_nts


class FeatureContainer(DataContainer):
    def __init__(self, name, bam_path, feature_path):
        print('Starting with features')
        super().__init__(name, bam_path)
        self.features = h5py.File(feature_path, 'r')

    def collect_nucleotides_aligned_to_mRNA_site(self, site, thresh_coverage=1, max_num_reads=1000, enforce_ref_5mer=False):
        site_nts = []
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
                        if enforce_ref_5mer and (query_5mer != site.ref_5mer):
                            continue

                        if query_name in self.features.keys():
                            site_nts.append(
                                Nucleotide(
                                    read_id=query_name,
                                    read_pos=query_position,
                                    ref_pos=site.start,
                                    pred_5mer=query_5mer,
                                    ref_5mer=site.ref_5mer,
                                    feature=np.array(self.features[query_name])[query_position],
                                    strand=site.strand
                                )
                            )
                            valid_counts+=1

                        if (max_num_reads > 0) and (len(site_nts) >= max_num_reads):
                            break
        if len(site_nts) > 0:
            self.nucleotides[site.ind] = site_nts