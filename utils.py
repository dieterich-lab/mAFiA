import os, argparse
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

class Args_Parser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--ref_file')
        self.add_argument('--max_num_reads', type=int, default=-1)
        self.add_argument('--min_coverage', type=int, default=0)
        self.add_argument('--enforce_ref_5mer', action='store_true')
        self.add_argument('--backbone_model_path')
        self.add_argument('--extraction_layer', default='convlayers.conv21')
        self.add_argument('--feature_width', type=int, default=0)
        self.add_argument('--classifier_type', default='logistic_regression')
        self.add_argument('--classifier_model_dir')

    def parse_and_print(self):
        self.args = self.parse_args()
        print('=========================================================')
        for k, v in vars(self.args).items():
            print('{} : {}'.format(k, v))
        print('=========================================================')

class Test_Args_Parser(Args_Parser):
    def __init__(self):
        super().__init__()
        self.add_argument('--test_bam_file')
        self.add_argument('--test_fast5_dir')
        self.add_argument('--outfile')

class Train_Args_Parser(Args_Parser):
    def __init__(self):
        super().__init__()
        self.add_argument('--unm_bam_file')
        self.add_argument('--unm_fast5_dir')
        self.add_argument('--mod_bam_file')
        self.add_argument('--mod_fast5_dir')
        self.add_argument('--scaler', default=None)

class mRNA_Test_Args_Parser(Test_Args_Parser):
    def __init__(self):
        super().__init__()
        self.add_argument('--mod_file')
        self.add_argument('--mod_prob_thresh', type=float, default=0.5)
        self.add_argument('--output_mod_probs', action='store_true')

class Output_Writer:
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

class mRNA_Output_Writer(Output_Writer):
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