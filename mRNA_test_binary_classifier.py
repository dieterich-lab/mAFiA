import os, sys, argparse
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import pandas as pd
from utils import load_reference, mRNA_output_writer
from data_containers import mRNA_site, mRNA_data_container
from feature_extractors import backbone_network
from feature_classifiers import load_motif_classifiers

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test_bam_file')
    parser.add_argument('--test_fast5_dir')
    parser.add_argument('--ref_file')
    parser.add_argument('--mod_file')
    parser.add_argument('--max_num_reads', type=int, default=-1)
    parser.add_argument('--min_coverage', type=int, default=0)
    parser.add_argument('--enforce_ref_5mer', action='store_true')
    parser.add_argument('--backbone_model_path')
    parser.add_argument('--extraction_layer', default='convlayers.conv21')
    parser.add_argument('--feature_width', type=int, default=0)
    parser.add_argument('--classifier_type', default='logistic_regression')
    parser.add_argument('--classifier_model_dir')
    parser.add_argument('--mod_prob_thresh', type=float, default=0.5)
    parser.add_argument('--outfile')
    parser.add_argument('--output_mod_probs', action='store_true')
    return parser.parse_args()

def main(args):
    test_container = mRNA_data_container('test', args.test_bam_file, args.test_fast5_dir)
    test_container.build_dict_read_ref()

    ivt_backbone = backbone_network(args.backbone_model_path, args.extraction_layer, args.feature_width)

    reference = load_reference(args.ref_file)

    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)

    writer = mRNA_output_writer(out_path=args.outfile, output_mod_probs=args.output_mod_probs)

    df_mod = pd.read_csv(args.mod_file)
    df_mod = df_mod.rename(columns={'Unnamed: 0': 'index'})
    for _, glori_row in df_mod[df_mod['index']>writer.last_ind].iterrows():
        this_mRNA_site = mRNA_site(glori_row, reference)
        if (this_mRNA_site.chr.isnumeric()==False) and (this_mRNA_site.chr not in ['X', 'Y']):
            continue
        if this_mRNA_site.ref_motif not in motif_classifiers.keys():
            continue

        test_container.collect_nucleotides_aligned_to_mRNA_site(
            ivt_backbone, this_mRNA_site,
            thresh_coverage=args.min_coverage,
            max_num_reads=args.max_num_reads,
            enforce_ref_5mer=args.enforce_ref_5mer
        )

        if len(test_container.nucleotides.get(this_mRNA_site.ind, [])) > args.min_coverage:
            print('=========================================================')
            this_mRNA_site.print()
            mod_ratio = motif_classifiers[this_mRNA_site.ref_motif].test(test_container.nucleotides[this_mRNA_site.ind])
            print('=========================================================\n')
            df_nts = test_container.flush_nts_to_dataframe()
            writer.update_df_out(glori_row, df_nts, mod_ratio)
            writer.write_df()
        else:
            test_container.nucleotides.clear()
    print('Total {} sites written to {}'.format(writer.site_counts, writer.out_path))

if __name__ == "__main__":
    main(parse_args())
