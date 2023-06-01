import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
from utils import test_args_parser, load_reference, parse_motif_dims, output_writer
from data_containers import oligo_data_container
from feature_extractors import backbone_network
from feature_classifiers import load_motif_classifiers

parser = test_args_parser()
parser.parse_and_print()

def main(args):
    test_container = oligo_data_container('test', args.test_bam_file, args.test_fast5_dir)
    test_container.build_dict_read_ref()

    ivt_backbone = backbone_network(args.backbone_model_path, args.extraction_layer, args.feature_width)

    test_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)

    reference = load_reference(args.ref_file)
    motif_dims = parse_motif_dims(reference)

    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)

    writer = output_writer(out_path=args.outfile)

    for motif_ind, motif, block_size, block_center in motif_dims:
        if motif not in motif_classifiers.keys():
            print('No classifier available for {}. Skipping...'.format(motif))
            continue

        test_container.collect_motif_nucleotides(motif_ind, motif, reference, block_size=block_size, block_center=block_center, enforce_ref_5mer=args.enforce_ref_5mer)

        if len(test_container.nucleotides[motif])<args.min_coverage:
            print('Insufficient coverage {}. Skipping...'.format(len(test_container.nucleotides[motif])))
            continue

        _ = motif_classifiers[motif].test(test_container.nucleotides[motif])

        df_nts = test_container.flush_nts_to_dataframe()
        writer.update_df_out(df_nts)
        writer.write_df()
    print('Total number of nucleotides tested {}'.format(len(writer.df_out)), flush=True)

if __name__ == "__main__":
    main(parser.args)
