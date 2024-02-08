import os
import time
import pandas as pd
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/mAFiA_dev'))
from mAFiA.arg_parsers import mRNATestArgsParser
from mAFiA.data_containers import MultiReadContainer
from mAFiA.feature_extractors import BackboneNetwork
from mAFiA.feature_classifiers import load_motif_classifiers
from mAFiA.output_writers import SAMWriter


def main():
    tic = time.time()

    parser = mRNATestArgsParser()
    parser.parse_and_print()
    args = parser.args

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    test_container = MultiReadContainer('test', args.bam_file, args.fast5_dir)
    test_container.build_dict_read_ref()
    ivt_backbone = BackboneNetwork(args.backbone_model_path, args.extraction_layer, args.feature_width, args.batchsize)
    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)
    sam_writer = SAMWriter(in_bam_path=args.bam_file, out_sam_path=os.path.join(args.out_dir, 'mAFiA.reads.sam'))
    df_mod = pd.read_csv(args.mod_file, sep='\t', dtype={'chrom': str, 'chromStart': int, 'chromEnd': int})
    df_mod_avail = df_mod[df_mod['ref5mer'].isin(motif_classifiers.keys())]

    test_container.process_reads(ivt_backbone, df_mod_avail, motif_classifiers, sam_writer)

    print(f'Total {sam_writer.read_counts} mod. reads written to {sam_writer.out_sam_path}')
    toc = time.time()
    print('Finished in {:.1f} mins'.format((toc - tic) / 60), flush=True)


if __name__ == "__main__":
    main()
