import os
import time
import pandas as pd
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/mAFiA_dev'))
from mAFiA.arg_parsers import mRNATestArgsParser
from mAFiA.data_containers import MultiReadContainer
from mAFiA.feature_extractors import BackboneNetwork
from mAFiA.feature_classifiers import load_multimod_motif_classifiers
from mAFiA.output_writers import SAMWriter
import pysam
from joblib import Parallel, delayed

def process_bam(in_bam_file, out_sam_file, args):
    test_container = MultiReadContainer('test', in_bam_file, args.fast5_dir)
    test_container.build_dict_read_ref()
    ivt_backbone = BackboneNetwork(args.backbone_model_path, args.extraction_layer, args.feature_width, args.batchsize)
    multimod_motif_classifiers = load_multimod_motif_classifiers(args.classifier_model_dir)
    sam_writer = SAMWriter(in_bam_path=in_bam_file, out_sam_path=out_sam_file)
    df_mod = pd.read_csv(args.mod_file, sep='\t', dtype={'chrom': str, 'chromStart': int, 'chromEnd': int})
    all_motifs = [motif for mod, motif_classifiers in multimod_motif_classifiers.items() for motif in motif_classifiers.keys()]
    df_mod_avail = df_mod[df_mod['ref5mer'].isin(all_motifs)]
    test_container.process_reads(ivt_backbone, df_mod_avail, multimod_motif_classifiers, sam_writer)

    # print(f'Total {sam_writer.read_counts} mod. reads written to {sam_writer.out_sam_path}')

    return sam_writer


def split_bam_file(in_bam_file, out_dir, num_jobs):
    out_bam_pattern = os.path.join(out_dir, 'temp%d.bam')

    # total_num_alignments = int(pysam.idxstats(in_bam_file).rstrip('\n').split('\t')[2])
    total_num_alignments = sum([int(l.split('\t')[2]) for l in pysam.idxstats(in_bam_file).split('\n') if len(l)])
    print(f'Splitting {total_num_alignments} reads into {num_jobs} chunks...', flush=True)
    chunk_size = total_num_alignments // num_jobs
    chunk = 0
    reads_in_this_chunk = 0
    old_name = None

    bam_files = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
        outfile = pysam.AlignmentFile(out_bam_pattern % chunk, 'wb', template=in_bam)
        bam_files.append(out_bam_pattern % chunk)
        for read in in_bam.fetch(until_eof=True):
            if old_name != read.query_name and reads_in_this_chunk > chunk_size:
                reads_in_this_chunk = 0
                chunk += 1
                outfile.close()
                outfile = pysam.AlignmentFile(out_bam_pattern % chunk, 'wb', template=in_bam)
                bam_files.append(out_bam_pattern % chunk)

            outfile.write(read)
            old_name = read.query_name
            reads_in_this_chunk += 1

    for this_bam_file in bam_files:
        pysam.index(this_bam_file)

    return bam_files


def main():
    tic = time.time()

    parser = mRNATestArgsParser()
    parser.parse_and_print()
    args = parser.args

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    split_bam_files = split_bam_file(args.bam_file, args.out_dir, args.num_jobs)
    out_sam_files = [in_bam_file.replace('.bam', 'mAFiA.reads.sam') for in_bam_file in split_bam_files]

    sam_writers = Parallel(n_jobs=args.num_jobs)(delayed(process_bam)(split_bam_files[i], out_sam_files[i], args) for i in range(len(split_bam_files)))
    total_num_reads_written = sum([this_writer.read_counts for this_writer in sam_writers])
    print(f'Total {total_num_reads_written} mod. reads written')

    pysam.merge('-f', '-o', os.path.join(args.out_dir, 'mAFiA.reads.bam'), *out_sam_files)
    pysam.index(os.path.join(args.out_dir, 'mAFiA.reads.bam'))

    toc = time.time()
    print('Finished in {:.1f} mins'.format((toc - tic) / 60), flush=True)


if __name__ == "__main__":
    main()
