import argparse

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
            print(f'{k} : {v}')
        print('=========================================================')

class Test_Args_Parser(Args_Parser):
    def __init__(self):
        super().__init__()
        self.add_argument('--test_bam_file')
        self.add_argument('--test_fast5_dir')
        self.add_argument('--out_dir')

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