import os
import pandas as pd
import pysam
import numpy as np

mAFiA_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
    'coverage',
    'modRatio',
    'confidence'
]

KS_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer',
    'ks_stat',
    'pval',
    'coverage_1',
    'coverage_2',
    'modRatio_1',
    'modRatio_2',
    'delta'
]

class DataframeWriter:
    def __init__(self, out_path, fmt_precision=6):
        self.out_path = out_path
        self.fmt_precision = fmt_precision
        outdir = os.path.dirname(self.out_path)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        # self.df_out = pd.DataFrame()

    def update_df_out(self, new_df):
        self.df_out = pd.concat([self.df_out, new_df]).round(self.fmt_precision)

    # def write_df(self, empty=False):
    #     self.df_out.to_csv(self.out_path, sep='\t', index=False)
    #     if empty:
    #         self.df_out = pd.DataFrame()

    def write_df(self, empty=False, bed_fields=mAFiA_fields):
        out_df = pd.DataFrame(self.out_rows, columns=bed_fields)
        if not os.path.exists(self.out_path):
            out_df.to_csv(self.out_path, sep='\t', index=False, header=True)
        else:
            out_df.to_csv(self.out_path, sep='\t', index=False, header=False, mode='a')
        if empty:
            self.out_rows = []


class SiteWriter(DataframeWriter):
    def __init__(self, out_path):
        super().__init__(out_path)

        # if os.path.exists(out_path):
        #     self.df_out = pd.read_csv(out_path, sep='\t')
        #     self.site_counts = len(self.df_out['index'].unique())
        #     self.site_counts = len(self.df_out)
        #     if self.site_counts > 0:
        #         self.last_ind = self.df_out.tail(1)['index'].values[0]
        #         print(f'Restarting from {out_path}, index {self.last_ind}, {self.site_counts} sites')
        #         return
        # else:
        # self.df_out = pd.DataFrame()
        self.out_rows = []
        self.site_counts = 0
        # self.last_ind = -1

    # def update_df_out(self, glori, nts, pred_ratio):
    #     df_glori = pd.concat([glori.to_frame().T] * len(nts), ignore_index=True)
    #     df_glori_nts = pd.concat([df_glori, nts], axis=1)
    #     df_glori_nts['frequency'] = round(pred_ratio*100.0)
    #     self.site_counts += 1
    #     self.df_out = pd.concat([self.df_out, df_glori_nts])

    def update_sites(self, in_row, cov, ratio, conf, ref_5mer, train_5mer=None):
        in_row['coverage'] = cov
        in_row['modRatio'] = round(ratio*100.0, 1)
        in_row['confidence'] = round(conf*100.0, 1)
        in_row['ref5mer'] = ref_5mer
        if train_5mer:
            in_row['train5mer'] = train_5mer
        # self.df_out = pd.concat([self.df_out, pd.DataFrame(in_row).T])
        self.out_rows.append(in_row)
        self.site_counts += 1


class KSWriter(DataframeWriter):
    def __init__(self, out_path):
        super().__init__(out_path)
        self.out_rows = []
        self.site_counts = 0

    def update_sites(self, in_row, ks_stat, pval, coverage_1, coverage_2):
        in_row['ks_stat'] = round(ks_stat, 4)
        in_row['pval'] = round(pval, 4)
        in_row['coverage_1'] = coverage_1
        in_row['coverage_2'] = coverage_2
        self.out_rows.append(in_row)
        self.site_counts += 1

    def write_df(self):
        super(KSWriter, self).write_df(empty=True, bed_fields=KS_fields)


class BAMWriter:
    def __init__(self, in_bam_path, out_bam_path):
        self.in_bam_path = in_bam_path
        self.out_bam_path = out_bam_path
        self.dict_read_mod = {}
        self.read_counts = 0

    def build_dict_read_mod(self, site_nts):
        all_nts = [nt for k, v in site_nts.items() for nt in v]
        for this_nt in all_nts:
            if this_nt.mod_prob>0:
                if this_nt.read_id not in self.dict_read_mod.keys():
                    self.dict_read_mod[this_nt.read_id] = []
                self.dict_read_mod[this_nt.read_id].append((this_nt.read_pos, this_nt.strand, this_nt.mod_prob, this_nt.pred_5mer, this_nt.ref_5mer))

    def generate_mm_ml_tags(self, read_mods, mod_base='N', mod_code='21891'):
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

    def write_bam_with_mm_ml_tags(self, container, mod_base, mod_code):
        self.build_dict_read_mod(container.nucleotides)
        with pysam.Samfile(self.in_bam_path, "rb") as fi:
            with pysam.Samfile(self.out_bam_path, "wb", template=fi) as fo:
                for this_read in fi.fetch():
                    this_read_mods = self.dict_read_mod.get(this_read.query_name)
                    if this_read_mods:
                        this_read_mods.sort(key=lambda x: x[0])
                        mm, ml = self.generate_mm_ml_tags(this_read_mods, mod_base, str(mod_code))
                        this_read.set_tag('MM', mm)
                        this_read.set_tag('ML', ml)
                    fo.write(this_read)
                    self.read_counts += 1


dict_mod_code = {
    'm6A': '21891',
    'psi': '17802',
    'Gm': '19229'
}


class SAMWriter:
    def __init__(self, in_bam_path, out_sam_path):
        self.in_bam_path = in_bam_path
        self.out_sam_path = out_sam_path
        self.fi = pysam.Samfile(in_bam_path, "rb")
        self.read_counts = 0

    def open(self):
        self.fo = pysam.Samfile(self.out_sam_path, "w", template=self.fi)

    def get_processed_reads(self):
        processed_reads = []
        with pysam.Samfile(self.out_sam_path, "r") as out_sam:
            try:
                for read in out_sam.fetch():
                    processed_reads.append(read)
            except:
                pass
        return processed_reads[::-1]

    def get_processed_read_ids(self):
        processed_read_ids = []
        with pysam.Samfile(self.out_sam_path, "r") as out_sam:
            try:
                for read in out_sam.fetch():
                    processed_read_ids.append(read.query_name)
            except:
                pass
        if len(processed_read_ids)>0:
            self.read_counts = len(processed_read_ids) - 1
            os.system('sed -i "$ d" {0}'.format(self.out_sam_path))
        return processed_read_ids[:-1]

    def build_dict_read_mod(self, read_nts):
        all_nts = [nt for k, v in read_nts.items() for nt in v]
        for this_nt in all_nts:
            if this_nt.mod_prob>0:
                if this_nt.read_id not in self.dict_read_mod.keys():
                    self.dict_read_mod[this_nt.read_id] = []
                self.dict_read_mod[this_nt.read_id].append((this_nt.read_pos, this_nt.strand, this_nt.mod_prob, this_nt.pred_5mer, this_nt.ref_5mer))

    def generate_mm_ml_tags(self, read_mods, mod_base, mod_code):
        dists = [read_mods[0][0]] + list(np.diff([mod[0] for mod in read_mods])-1)

        unique_strands = np.unique([mod[1] for mod in read_mods])
        mod_probs = [mod[2] for mod in read_mods]
        rescaled_mod_probs = [round(mp*255.0) for mp in mod_probs]
        if len(unique_strands)==1:
            strand = unique_strands[0]
        else:
            print('Warning: mixed strands!!!')
        mm_tag = mod_base + strand + mod_code + ',' + ','.join([str(d) for d in dists]) + ';'
        ml_tag = rescaled_mod_probs
        return mm_tag, ml_tag

    def write_read(self, read, read_nts, mod_base = 'N'):
        full_mm = ''
        full_ml = []
        for this_mod, this_mod_nts in read_nts.items():
            read_mods = [
                (this_nt.read_pos, this_nt.strand, this_nt.mod_prob, this_nt.pred_5mer, this_nt.ref_5mer)
                for this_nt in this_mod_nts]
            this_mod_code = dict_mod_code[this_mod]
            if read_mods:
                read_mods.sort(key=lambda x: x[0])
                this_mod_mm, this_mod_ml = self.generate_mm_ml_tags(read_mods, mod_base, str(this_mod_code))
                full_mm += this_mod_mm
                full_ml += this_mod_ml
        if len(full_mm)>0:
            read.set_tag('MM', full_mm)
            read.set_tag('ML', full_ml)
        self.fo.write(read)
        self.read_counts += 1

    def write_reads(self, in_reads_mod_nts):
        for write_read, write_mod_nts in in_reads_mod_nts:
            self.write_read(write_read, write_mod_nts)

    def close(self):
        self.fi.close()
        self.fo.close()