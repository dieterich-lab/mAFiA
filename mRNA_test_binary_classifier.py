import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import pandas as pd
from Bio import SeqIO
from arg_parsers import mRNA_Test_Args_Parser, mRNA_Output_Writer
from data_containers import mRNA_Site, mRNA_Data_Container
from feature_extractors import Backbone_Network
from feature_classifiers import load_motif_classifiers

parser = mRNA_Test_Args_Parser()
parser.parse_and_print()

def load_genome_reference(ref_file):
    print('Parsing genome reference {}...'.format(os.path.basename(ref_file)))
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
            ref[record.id] = str(record.seq)
        elif record.id=='MT':
            ref['M'] = str(record.seq)
    return ref

def main(args):
    test_container = mRNA_Data_Container('test', args.test_bam_file, args.test_fast5_dir)
    test_container.build_dict_read_ref()

    ivt_backbone = Backbone_Network(args.backbone_model_path, args.extraction_layer, args.feature_width)
    reference = load_genome_reference(args.ref_file)
    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)
    writer = mRNA_Output_Writer(out_path=args.outfile, output_mod_probs=args.output_mod_probs)

    df_mod = pd.read_csv(args.mod_file)
    df_mod = df_mod.rename(columns={'Unnamed: 0': 'index'})
    for _, glori_row in df_mod[df_mod['index']>writer.last_ind].iterrows():
        this_mRNA_site = mRNA_Site(glori_row, reference)
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
