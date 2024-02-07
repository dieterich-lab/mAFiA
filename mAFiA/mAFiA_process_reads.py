import os
import time
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/mAFiA_dev'))
from mAFiA.arg_parsers import mRNATestArgsParser
from mAFiA.data_containers import mRNASite, MultiReadContainer
from mAFiA.feature_extractors import BackboneNetwork
from mAFiA.feature_classifiers import load_motif_classifiers
from mAFiA.output_writers import SAMWriter


def load_genome_reference(ref_file):
    print(f'Parsing genome reference {os.path.basename(ref_file)}...')
    ref = {}
    # for record in SeqIO.parse(ref_file, 'fasta'):
    #     if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
    #         ref[record.id] = str(record.seq)
    #     elif record.id=='MT':
    #         ref['M'] = str(record.seq)
    for record in SeqIO.parse(ref_file, 'fasta'):
        ref[record.id] = str(record.seq)

    return ref


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

    test_container.process_reads(ivt_backbone, df_mod, motif_classifiers, sam_writer)



# for _, row in tqdm(list(df_mod.iterrows())):
#     this_mRNA_site = mRNASite(row, reference)
#     if this_mRNA_site.ref_5mer not in motif_classifiers.keys():
#         continue
#
#     if args.features_file is not None:
#         test_container.collect_nucleotides_aligned_to_mRNA_site(
#             site=this_mRNA_site,
#             thresh_coverage=args.min_coverage,
#             max_num_reads=args.max_num_reads,
#             enforce_ref_5mer=args.enforce_ref_5mer
#         )
#     else:
#         test_container.collect_nucleotides_aligned_to_mRNA_site(
#             ivt_backbone, this_mRNA_site,
#             thresh_coverage=args.min_coverage,
#             max_num_reads=args.max_num_reads,
#             enforce_ref_5mer=args.enforce_ref_5mer
#         )
#
#     this_site_coverage = len(test_container.nucleotides.get(this_mRNA_site.ind, []))
#     if this_site_coverage >= args.min_coverage:
#         print('=========================================================')
#         this_mRNA_site.print()
#         this_site_mod_ratio = motif_classifiers[this_mRNA_site.ref_5mer].test(
#             test_container.nucleotides[this_mRNA_site.ind]
#         )
#         print('=========================================================\n')
#         site_writer.update_site_df(row, this_site_coverage, this_site_mod_ratio, this_mRNA_site.ref_5mer)
#         site_writer.write_df()
# print(f'Total {site_writer.site_counts} mod. sites written to {site_writer.out_path}')
#
# bam_writer.write_bam_with_mm_ml_tags(test_container)
# print(f'Total {bam_writer.read_counts} mod. reads written to {bam_writer.out_bam_path}')
# toc = time.time()
# print('Finished in {:.1f} mins'.format((toc - tic) / 60), flush=True)


if __name__ == "__main__":
    main()
