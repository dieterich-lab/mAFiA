import os
HOME = os.path.expanduser('~')
import sys
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
import argparse
import pandas as pd
from Bio.Seq import Seq
from utils import load_reference
from data_containers import mRNA_site, mRNA_data_container
from feature_extractors import backbone_network
from feature_classifiers import load_motif_classifiers

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
args = parser.parse_args()

outdir = os.path.dirname(args.outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

test_container = mRNA_data_container('test', args.test_bam_file, args.test_fast5_dir)
test_container.build_dict_read_ref()

ivt_backbone = backbone_network(args.backbone_model_path, args.extraction_layer, args.feature_width)

reference = load_reference(args.ref_file)

motif_classifiers = load_motif_classifiers(args.classifier_model_dir)

### check for existing output ###
# restart = False
# if os.path.exists(args.outfile):
#     df_out = pd.read_csv(args.outfile, sep='\t')
#     site_counts = len(df_out)
#     if site_counts>0:
#         print('Restarting from {} with {} counts'.format(args.outfile, site_counts))
#         last_row = df_out.tail(1)
#         last_chr = last_row['Chr'].values[0].lstrip('chr')
#         last_start = last_row['Sites'].values[0] - 1  # 0-based
#         restart = True
# else:
#     df_out = pd.DataFrame()
#     site_counts = 0
#     print('Starting from scratch')
#     restart = False

df_mod = pd.read_csv(args.mod_file)
df_mod = df_mod.rename(columns={'Unnamed: 0': 'index'})
df_out = pd.DataFrame()
site_counts = 0
for _, glori_row in df_mod.iterrows():
    this_mRNA_site = mRNA_site(glori_row, reference)
    if (this_mRNA_site.chr.isnumeric()==False) and (this_mRNA_site.chr not in ['X', 'Y']):
        continue
    if this_mRNA_site.ref_motif not in motif_classifiers.keys():
        continue

    # if restart:
    #     if (chr==last_chr) and (start==last_start):
    #         print('Last row in file: chr{}, pos{}'.format(chr, start))
    #         restart = False
    #     else:
    #         print('Skipping chr{}, pos{}'.format(chr, start))
    #     continue

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
        df_nts = test_container.flush_nts_to_dataframe()
        print('=========================================================\n')

        df_glori_row = pd.concat([glori_row.to_frame().T] * len(df_nts), ignore_index=True)
        df_glori_nts = pd.concat([df_glori_row, df_nts], axis=1)
        df_glori_nts['pred_mod_ratio'] = round(mod_ratio, 3)
        df_out = pd.concat([df_out, df_glori_nts])
        df_out.to_csv(args.outfile, sep='\t', index=False)
        site_counts += 1
    else:
        _ = test_container.nucleotides.clear()
print('Total {} sites written to {}'.format(site_counts, args.outfile))