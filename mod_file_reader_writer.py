import os
import pandas as pd

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

        # if os.path.exists(out_path):
        #     self.df_out = pd.read_csv(out_path, sep='\t')
        #     self.site_counts = len(self.df_out['index'].unique())
        #     if self.site_counts > 0:
        #         self.last_ind = self.df_out.tail(1)['index'].values[0]
        #         print('Restarting from {}, index {}, {} sites'.format(out_path, self.last_ind, self.site_counts))
        #         return
        self.df_out = pd.DataFrame()
        self.site_counts = 0
        # self.last_ind = -1
        print('Starting from scratch')

    def update_df_out(self, glori, nts, pred_ratio):
        df_glori = pd.concat([glori.to_frame().T] * len(nts), ignore_index=True)
        df_glori_nts = pd.concat([df_glori, nts], axis=1)
        df_glori_nts['frequency'] = round(pred_ratio*100.0)
        self.site_counts += 1
        self.df_out = pd.concat([self.df_out, df_glori_nts])

    def update_site_df(self, in_row, cov, ratio, ref_5mer):
        in_row['coverage'] = cov
        in_row['frequency'] = round(ratio*100.0)
        in_row['ref5mer'] = ref_5mer
        self.site_counts += 1
        self.df_out = pd.concat([self.df_out, pd.DataFrame(in_row).T])
