import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from arg_parsers import mRNA_Test_Args_Parser
from data_containers import mRNA_Site, mRNA_Data_Container
from feature_extractors import Backbone_Network
from feature_classifiers import load_motif_classifiers
from output_writers import Site_Writer, BAM_Writer

parser = mRNA_Test_Args_Parser()
parser.parse_and_print()

def load_genome_reference(ref_file):
    print(f'Parsing genome reference {os.path.basename(ref_file)}...')
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        if (record.id.isnumeric()) or (record.id in ['X', 'Y']):
            ref[record.id] = str(record.seq)
        elif record.id=='MT':
            ref['M'] = str(record.seq)
    return ref

def main(args):
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    test_container = mRNA_Data_Container('test', args.test_bam_file, args.test_fast5_dir)
    test_container.build_dict_read_ref()

    ivt_backbone = Backbone_Network(args.backbone_model_path, args.extraction_layer, args.feature_width)
    reference = load_genome_reference(args.ref_file)
    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)
    site_writer = Site_Writer(out_path=os.path.join(args.out_dir, 'site_mods.bed'))
    bam_writer = BAM_Writer(in_bam_path=args.test_bam_file, out_bam_path=os.path.join(args.out_dir, 'read_mods.bam'))

    df_mod = pd.read_csv(args.mod_file, sep='\t')
    for _, row in tqdm(list(df_mod.iterrows())):
        this_mRNA_site = mRNA_Site(row, reference)
        if this_mRNA_site.ref_5mer not in motif_classifiers.keys():
            continue

        test_container.collect_nucleotides_aligned_to_mRNA_site(
            ivt_backbone, this_mRNA_site,
            thresh_coverage=args.min_coverage,
            max_num_reads=args.max_num_reads,
            enforce_ref_5mer=args.enforce_ref_5mer
        )

        this_site_coverage = len(test_container.nucleotides.get(this_mRNA_site.ind, []))
        if this_site_coverage > args.min_coverage:
            print('=========================================================')
            this_mRNA_site.print()
            this_site_mod_ratio = motif_classifiers[this_mRNA_site.ref_5mer].test(test_container.nucleotides[this_mRNA_site.ind])
            print('=========================================================\n')
            # df_nts = test_container.flush_nts_to_dataframe()
            # df_nts = test_container.flush_nts_to_dataframe()
            site_writer.update_site_df(row, this_site_coverage, this_site_mod_ratio, this_mRNA_site.ref_5mer)
            site_writer.write_df()
        # else:
            # test_container.nucleotides.clear()
    print(f'Total {site_writer.site_counts} mod. sites written to {site_writer.out_path}')

    bam_writer.write_bam_with_mm_ml_tags(test_container)
    print(f'Total {bam_writer.read_counts} mod. reads written to {bam_writer.out_bam_path}')

if __name__ == "__main__":
    main(parser.args)
