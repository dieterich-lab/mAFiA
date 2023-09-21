import os
import pandas as pd
import pysam
import numpy as np

class DataframeWriter:
    def __init__(self, out_path, fmt_precision=6):
        self.out_path = out_path
        self.fmt_precision = fmt_precision
        outdir = os.path.dirname(self.out_path)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        self.df_out = pd.DataFrame()

    def update_df_out(self, new_df):
        self.df_out = pd.concat([self.df_out, new_df]).round(self.fmt_precision)

    def write_df(self):
        self.df_out.to_csv(self.out_path, sep='\t', index=False)

class SiteWriter(DataframeWriter):
    def __init__(self, out_path):
        super().__init__(out_path)

        # if os.path.exists(out_path):
        #     self.df_out = pd.read_csv(out_path, sep='\t')
        #     self.site_counts = len(self.df_out['index'].unique())
        #     if self.site_counts > 0:
        #         self.last_ind = self.df_out.tail(1)['index'].values[0]
        #         print(f'Restarting from {out_path}, index {self.last_ind}, {self.site_counts} sites')
        #         return
        self.df_out = pd.DataFrame()
        self.site_counts = 0
        # self.last_ind = -1

    def update_df_out(self, glori, nts, pred_ratio):
        df_glori = pd.concat([glori.to_frame().T] * len(nts), ignore_index=True)
        df_glori_nts = pd.concat([df_glori, nts], axis=1)
        df_glori_nts['frequency'] = round(pred_ratio*100.0)
        self.site_counts += 1
        self.df_out = pd.concat([self.df_out, df_glori_nts])

    def update_site_df(self, in_row, cov, ratio, ref_5mer, train_5mer=None):
        in_row['coverage'] = cov
        in_row['modRatio'] = round(ratio*100.0)
        in_row['ref5mer'] = ref_5mer
        if train_5mer:
            in_row['train5mer'] = train_5mer
        self.site_counts += 1
        self.df_out = pd.concat([self.df_out, pd.DataFrame(in_row).T])

class BAMWriter:
    def __init__(self, in_bam_path, out_bam_path):
        self.in_bam_path = in_bam_path
        self.out_bam_path = out_bam_path
        self.dict_read_mod = {}
        self.read_counts = 0
        self.dict_ChEBI = {
            'm6A': '21891',
            'Gm': '19229'
        }

    def build_dict_read_mod(self, site_nts):
        all_nts = [nt for k, v in site_nts.items() for nt in v]
        for this_nt in all_nts:
            if this_nt.mod_prob>0:
                if this_nt.read_id not in self.dict_read_mod.keys():
                    self.dict_read_mod[this_nt.read_id] = []
                self.dict_read_mod[this_nt.read_id].append((this_nt.read_pos, this_nt.strand, this_nt.mod_prob, this_nt.pred_5mer, this_nt.ref_5mer))

    def generate_mm_ml_tags(self, read_mods, mod_base='N', mod_name='m6A'):
        mod_code = self.dict_ChEBI[mod_name]
        dists = [read_mods[0][0]] + list(np.diff([mod[0] for mod in read_mods])-1)
        unique_strands = np.unique([mod[1] for mod in read_mods])
        mod_probs = [mod[2] for mod in read_mods]
        rescaled_mod_probs = [round(mp*255.0) for mp in mod_probs]
        if len(unique_strands)==1:
            strand = unique_strands[0]
        else:
            print('Warning: mixed strands!!!')
        # mm_tag = 'MM:Z:N' + strand + str(mod_code) + ',' + ','.join([str(d) for d in dists]) + ';'
        # ml_tag = 'ML:B:C,' + ','.join([str(p) for p in rescaled_mod_probs])
        mm_tag = mod_base + strand + mod_code + ',' + ','.join([str(d) for d in dists]) + ';'
        ml_tag = rescaled_mod_probs
        return mm_tag, ml_tag

    def write_bam_with_mm_ml_tags(self, container, site):
        self.build_dict_read_mod(container.nucleotides)
        with pysam.Samfile(self.in_bam_path, "rb") as fi:
            with pysam.Samfile(self.out_bam_path, "wb", template=fi) as fo:
                for this_read in fi.fetch():
                    this_read_mods = self.dict_read_mod.get(this_read.query_name)
                    if this_read_mods:
                        this_read_mods.sort(key=lambda x: x[0])
                        mm, ml = self.generate_mm_ml_tags(this_read_mods, site.mod_name)
                        this_read.set_tag('MM', mm)
                        this_read.set_tag('ML', ml)
                        fo.write(this_read)
                        self.read_counts += 1
