import argparse


class ArgsParser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--ref_file')
        self.add_argument('--max_num_reads', type=int, default=-1)
        self.add_argument('--min_coverage', type=int, default=1)
        self.add_argument('--enforce_ref_5mer', action='store_true')
        self.add_argument('--backbone_model_path')
        self.add_argument('--extraction_layer', default='convlayers.conv21')
        self.add_argument('--feature_width', type=int, default=0)
        self.add_argument('--classifier_type', default='logistic_regression')
        self.add_argument('--classifier_model_dir')
        self.args = None

    def parse_and_print(self):
        self.args = self.parse_args()
        print('\n=========================================================')
        for k, v in vars(self.args).items():
            print(f'{k} : {v}')
        print('=========================================================')


class TestArgsParser(ArgsParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--bam_file')
        self.add_argument('--fast5_dir')
        self.add_argument('--out_dir')
        self.add_argument('--out_filename', default=None)
        self.add_argument('--batchsize', type=int, default=1024)

class TrainArgsParser(ArgsParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--unm_bam_file')
        self.add_argument('--unm_fast5_dir')
        self.add_argument('--mod_bam_file')
        self.add_argument('--mod_fast5_dir')
        self.add_argument('--scaler', default=None)

class mRNATestArgsParser(TestArgsParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--features_file', default=None)
        self.add_argument('--mod_file')
        self.add_argument('--mod_prob_thresh', type=float, default=0.5)
        self.add_argument('--output_mod_probs', action='store_true')
        self.add_argument('--chunk_size', type=int, default=10000)
        self.add_argument('--num_jobs', type=int, default=16)

class KSTestArgsParser(ArgsParser):
    def __init__(self):
        super().__init__()
        self.add_argument('--bam_file_1')
        self.add_argument('--bam_file_2')
        self.add_argument('--bed_file_1')
        self.add_argument('--bed_file_2')
        self.add_argument('--min_confidence', type=float, default=0.0)
        self.add_argument('--mod_file')
        self.add_argument('--out_dir')
        self.add_argument('--out_filename', default=None)
        self.add_argument('--chunk_size', type=int, default=1000)
        self.add_argument('--num_jobs', type=int, default=2)

