import os
import time
import pandas as pd
import pysam
from statistics import mean
from tqdm import tqdm
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/mAFiA_dev'))
from mAFiA.arg_parsers import mRNATestArgsParser
from mAFiA.output_writers import SiteWriter


dict_mod_code = {
    'm6A': '21891',
    'psi': '17802'
}


def main():
    tic = time.time()

    parser = mRNATestArgsParser()
    parser.parse_and_print()
    args = parser.args

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    site_writer = SiteWriter(out_path=os.path.join(args.out_dir, 'mAFiA.sites.bed'))
    df_mod = pd.read_csv(args.mod_file, sep='\t', dtype={'chrom': str, 'chromStart': int, 'chromEnd': int})

    with pysam.Samfile(args.bam_file, 'r') as bam:
        for _, row in tqdm(list(df_mod.iterrows())):
            chrom = row['chrom']
            chromStart = row['chromStart']
            strand = row['strand']
            mod_type = row['name']
            mod_code = int(dict_mod_code[mod_type])
            flag_require = 0 if strand =='+' else 16
            for pileupcolumn in bam.pileup(chrom, chromStart, chromStart + 1, truncate=True, flag_require=flag_require):
                if pileupcolumn.reference_pos == chromStart:
                    this_site_coverage = pileupcolumn.get_num_aligned()
                    if this_site_coverage >= args.min_coverage:
                        mod_counts = []
                        for pileupread in pileupcolumn.pileups:
                            flag = pileupread.alignment.flag
                            query_position = pileupread.query_position
                            if query_position is None:
                                continue
                            if flag == 16:
                                query_position = pileupread.alignment.query_length - query_position - 1
                            mod_key = ('N', 0, mod_code) if flag==0 else ('N', 1, mod_code)
                            sel_tup = [tup for tup in pileupread.alignment.modified_bases_forward.get(mod_key, []) if tup[0]==query_position]
                            if len(sel_tup)==1:
                                mod_counts.append((sel_tup[0][1] / 255.0) >= args.mod_prob_thresh)
                        if len(mod_counts)>=args.min_coverage:
                            row['score'] = '.'
                            site_writer.update_sites(row, cov=len(mod_counts), ratio=mean(mod_counts), ref_5mer=row['ref5mer'])
            if site_writer.site_counts%args.chunk_size==0:
                site_writer.write_df(empty=True)
    site_writer.write_df(empty=True)
    print(f'Total {site_writer.site_counts} mod. sites written to {site_writer.out_path}')

    toc = time.time()
    print('Finished in {:.1f} mins'.format((toc - tic) / 60), flush=True)


if __name__ == "__main__":
    main()
