import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from arg_parsers import mRNATestArgsParser
from data_containers import mRNASite, FeatureContainer
from feature_classifiers import load_motif_classifiers
from output_writers import SiteWriter, BAMWriter

parser = mRNATestArgsParser()
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

    test_container = FeatureContainer('test', args.in_bam_file, args.in_feat_file)

    reference = load_genome_reference(args.ref_file)
    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)

    site_writer = SiteWriter(out_path=os.path.join(args.out_dir, 'mAFiA.sites.bed'))
    bam_writer = BAMWriter(in_bam_path=args.in_bam_file, out_bam_path=os.path.join(args.out_dir, 'mAFiA.reads.bam'))

    df_mod = pd.read_csv(args.mod_file, sep='\t')
    for _, row in tqdm(list(df_mod.iterrows())):
        this_mRNA_site = mRNASite(row, reference)

        test_container.collect_nucleotides_aligned_to_mRNA_site(
            site=this_mRNA_site,
            thresh_coverage=args.min_coverage,
            max_num_reads=args.max_num_reads,
            enforce_ref_5mer=args.enforce_ref_5mer
        )

        this_site_coverage = len(test_container.nucleotides.get(this_mRNA_site.ind, []))
        if this_site_coverage > args.min_coverage:
            for train_motif in motif_classifiers.keys():
                print('\n=========================================================')
                this_mRNA_site.print()
                print(f'Train motif {train_motif}')
                this_site_mod_ratio = motif_classifiers[train_motif].test(test_container.nucleotides[this_mRNA_site.ind])
                site_writer.update_sites(row, this_site_coverage, this_site_mod_ratio, this_mRNA_site.ref_5mer, train_motif)
                site_writer.write_df()
                print('=========================================================\n')
    print(f'Total {site_writer.site_counts} mod. sites written to {site_writer.out_path}')

    # bam_writer.write_bam_with_mm_ml_tags(test_container)
    # print(f'Total {bam_writer.read_counts} mod. reads written to {bam_writer.out_bam_path}')

if __name__ == "__main__":
    main(parser.args)
