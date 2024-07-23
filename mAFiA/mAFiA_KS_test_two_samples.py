import os
import time
import pandas as pd
import pysam
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/mAFiA_dev'))
from mAFiA.arg_parsers import KSTestArgsParser
from mAFiA.output_writers import KSWriter
from joblib import Parallel, delayed
import numpy as np
from scipy.stats import kstest


dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
    'Gm': 19229
}


def collect_mod_prob(in_bam_file, in_site):
    in_chrom, in_chromStart, in_chromEnd, in_name, in_strand = in_site[
        ['chrom', 'chromStart', 'chromEnd', 'name', 'strand']]
    out_mod_probs = []
    flag_require = 0 if in_strand == '+' else 16
    with pysam.Samfile(in_bam_file, 'rb') as in_bam:
        for pileupcolumn in in_bam.pileup(in_chrom, in_chromStart, in_chromEnd, truncate=True,
                                          flag_require=flag_require):
            if pileupcolumn.pos == in_chromStart:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        this_read_mod_probs = [score for (pos, score) in pileupread.alignment.modified_bases.get(
                            ('N', 0, dict_mod_code[in_name]), [])
                                               if pos == pileupread.query_position]
                        out_mod_probs.extend(this_read_mod_probs)
    return np.array(out_mod_probs) / 255.0


def KSTest_on_single_site(in_row, args):
    mod_probs_1 = collect_mod_prob(args.bam_file_1, in_row)
    mod_probs_2 = collect_mod_prob(args.bam_file_2, in_row)
    if (len(mod_probs_1)>=args.min_coverage) and (len(mod_probs_2)>=args.min_coverage):
        ks_stat, pval = kstest(mod_probs_1, mod_probs_2)
        return {'in_row': in_row, 'coverage_1': len(mod_probs_1), 'coverage_2': len(mod_probs_2),'ks_stat': ks_stat, 'pval': pval}
    else:
        return {}


def main():
    tic = time.time()

    parser = KSTestArgsParser()
    parser.parse_and_print()
    args = parser.args

    os.makedirs(args.out_dir, exist_ok=True)
    if args.out_filename is not None:
        site_writer = KSWriter(out_path=os.path.join(args.out_dir, args.out_filename))
    else:
        site_writer = KSWriter(out_path=os.path.join(args.out_dir, 'mAFiA.KSTest.bed'))
    df_mod = pd.read_csv(args.mod_file, sep='\t', dtype={'chrom': str, 'chromStart': int, 'chromEnd': int}, iterator=True, chunksize=args.chunk_size)

    for chunk in df_mod:
        sites = Parallel(n_jobs=args.num_jobs)(delayed(KSTest_on_single_site)(chunk.iloc[i], args) for i in range(len(chunk)))
        for this_site in sites:
            if this_site:
                site_writer.update_sites(**this_site)
        site_writer.write_df()
        print(f'{chunk.index[-1]+1} rows processed', flush=True)

    print(f'Total {site_writer.site_counts} mod. sites written to {site_writer.out_path}', flush=True)
    toc = time.time()
    print('Finished in {:.1f} mins'.format((toc - tic) / 60), flush=True)


if __name__ == "__main__":
    main()
