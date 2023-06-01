import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from ont_fast5_api.fast5_interface import get_fast5_file
from Bio import SeqIO

def med_mad(x, factor=1.4826):
    med = np.median(x)
    mad = np.median(np.absolute(x - med)) * factor
    return med, mad

def segment(seg, s):
    seg = np.concatenate((seg, np.zeros((-len(seg)%s))))
    nrows=((seg.size-s)//s)+1
    n=seg.strides[0]
    return np.lib.stride_tricks.as_strided(seg, shape=(nrows,s), strides=(s*n, n))

def get_norm_signal_from_read_id(id, index_paths):
    filepath = index_paths[id]
    f5 = get_fast5_file(filepath, mode="r")
    read = f5.get_read(id)
    signal = read.get_raw_data(scale=True)
    signal_start = 0
    signal_end = len(signal)
    med, mad = med_mad(signal[signal_start:signal_end])
    return (signal[signal_start:signal_end] - med) / mad

def index_fast5_files(f5_paths, bam=None):
    index_read_ids = {}
    query_names = []
    if bam is not None:
        query_names = [alignment.query_name for alignment in bam.fetch()]
    for f5_filepath in tqdm(f5_paths):
        try:
            f5 = get_fast5_file(f5_filepath, mode="r")
        except:
            print('Error reading {}!'.format(f5_filepath))
        else:
            read_ids = f5.get_read_ids()
            if len(query_names)>0:
                for read_id in read_ids:
                    if read_id in query_names:
                        index_read_ids[read_id] = f5_filepath
            else:
                for read_id in read_ids:
                    index_read_ids[read_id] = f5_filepath
    return index_read_ids


def load_reference(ref_file):
    print('Parsing reference {}...'.format(os.path.basename(ref_file)))
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        ref[record.id] = str(record.seq)

    return ref

def parse_motif_dims(ref):
    block_index_motif_size_center = []
    for k in ref.keys():
        if k.lstrip('block').split('_')[0] == '1':
            block_index = k.split('_')[1]
            block_seq = ref[k]
            block_size = len(block_seq)
            block_center = block_size // 2
            motif = block_seq[block_center - 2:block_center + 3]
            block_index_motif_size_center.append((block_index, motif, block_size, block_center))
    return block_index_motif_size_center

class output_writer:
    def __init__(self, out_path, output_mod_probs=True, fmt_precision=6):
        self.out_path = out_path
        self.output_mod_probs = output_mod_probs
        self.fmt_precision = fmt_precision
        outdir = os.path.dirname(self.out_path)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        self.df_out = pd.DataFrame()

    def update_df_out(self, nts):
        self.df_out = pd.concat([self.df_out, nts]).round(self.fmt_precision)

    def write_df(self):
        self.df_out.to_csv(self.out_path, sep='\t', index=False)

class mRNA_output_writer(output_writer):
    def __init__(self, out_path, output_mod_probs=True):
        super().__init__(out_path, output_mod_probs)

        if os.path.exists(out_path):
            self.df_out = pd.read_csv(out_path, sep='\t')
            self.site_counts = len(self.df_out['index'].unique())
            if self.site_counts > 0:
                self.last_ind = self.df_out.tail(1)['index'].values[0]
                print('Restarting from {}, index {}, {} sites'.format(out_path, self.last_ind, self.site_counts))
                return
        self.df_out = pd.DataFrame()
        self.site_counts = 0
        self.last_ind = -1
        print('Starting from scratch')

    def update_df_out(self, glori, nts, pred_ratio):
        df_glori = pd.concat([glori.to_frame().T] * len(nts), ignore_index=True)
        df_glori_nts = pd.concat([df_glori, nts], axis=1)
        df_glori_nts['pred_mod_ratio'] = pred_ratio
        self.site_counts += 1
        self.df_out = pd.concat([self.df_out, df_glori_nts])