import os, sys, argparse
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
from utils import load_reference, parse_motif_dims
from data_containers import oligo_data_container
from feature_extractors import backbone_network
from feature_classifiers import motif_classifier

parser = argparse.ArgumentParser()
parser.add_argument('--unm_bam_file')
parser.add_argument('--unm_fast5_dir')
parser.add_argument('--mod_bam_file')
parser.add_argument('--mod_fast5_dir')
parser.add_argument('--ref_file')
parser.add_argument('--max_num_reads', type=int, default=-1)
parser.add_argument('--min_coverage', type=int, default=1)
parser.add_argument('--enforce_ref_5mer', action='store_true')
parser.add_argument('--backbone_model_path')
parser.add_argument('--extraction_layer', default='convlayers.conv21')
parser.add_argument('--feature_width', type=int, default=0)
parser.add_argument('--scaler', default=None)
parser.add_argument('--classifier_type', default='logistic_regression')
parser.add_argument('--classifier_model_dir')
args = parser.parse_args()
if args.enforce_ref_5mer:
    classifier_model_dir = args.classifier_model_dir + '_enforceMotif'
else:
    classifier_model_dir = args.classifier_model_dir
os.makedirs(classifier_model_dir, exist_ok=True)

unm_container = oligo_data_container('unm', args.unm_bam_file, args.unm_fast5_dir)
mod_container = oligo_data_container('mod', args.mod_bam_file, args.mod_fast5_dir)

ivt_backbone = backbone_network(args.backbone_model_path, args.extraction_layer, args.feature_width)

unm_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)
mod_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)

reference = load_reference(args.ref_file)
motif_dims = parse_motif_dims(reference)

for motif_ind, motif, block_size, block_center in motif_dims:
    unm_container.collect_motif_nucleotides(motif_ind, motif, reference, block_size=block_size, block_center=block_center, enforce_ref_5mer=args.enforce_ref_5mer)
    mod_container.collect_motif_nucleotides(motif_ind, motif, reference, block_size=block_size, block_center=block_center, enforce_ref_5mer=args.enforce_ref_5mer)

    if min(len(unm_container.nucleotides[motif]), len(mod_container.nucleotides[motif]))>=args.min_coverage:
        this_motif_classifier = motif_classifier(motif=motif, classifier_type=args.classifier_type, scaler=args.scaler)
        this_motif_classifier.train(unm_container.nucleotides[motif], mod_container.nucleotides[motif])
        this_motif_classifier.save(
            out_model_path=os.path.join(classifier_model_dir, '{}.pkl'.format(motif)),
            draw_prc=True
        )
