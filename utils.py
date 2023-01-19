import numpy as np
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