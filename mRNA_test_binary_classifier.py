import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import pandas as pd
from utils import mRNA_test_args_parser, load_reference, mRNA_output_writer
from data_containers import mRNA_site, mRNA_data_container
from feature_extractors import backbone_network
from feature_classifiers import load_motif_classifiers

parser = mRNA_test_args_parser()
parser.parse_and_print()

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
    main(parser.args)
