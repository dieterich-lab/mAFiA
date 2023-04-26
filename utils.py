import numpy as np
from tqdm import tqdm
from ont_fast5_api.fast5_interface import get_fast5_file

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

def index_fast5_files(f5_paths, bam):
    index_read_ids = {}
    query_names = [alignment.query_name for alignment in bam.fetch()]
    for f5_filepath in tqdm(f5_paths):
        # if 'fail' in f5_filepath:
        #     print('Skipping {}'.format(f5_filepath))
        #     continue
        try:
            f5 = get_fast5_file(f5_filepath, mode="r")
            read_ids = f5.get_read_ids()
        except:
            print('Error reading {}!'.format(f5_filepath))
        else:
            for read_id in read_ids:
                if read_id in query_names:
                    index_read_ids[read_id] = f5_filepath
    return index_read_ids