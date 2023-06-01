import os
import numpy as np
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