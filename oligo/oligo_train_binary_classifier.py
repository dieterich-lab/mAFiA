import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
from mAFiA.arg_parsers import TrainArgsParser
from oligo_processors import OligoReferenceGenerator
from mAFiA.data_containers import OligoDataContainer
from mAFiA.feature_extractors import BackboneNetwork
from mAFiA.feature_classifiers import MotifClassifier

parser = TrainArgsParser()
parser.parse_and_print()

### todo: collect unm and mod features in same BAM file ###

def main(args):
    os.makedirs(args.classifier_model_dir, exist_ok=True)

    unm_container = OligoDataContainer('unm', args.unm_bam_file, args.unm_fast5_dir)
    mod_container = OligoDataContainer('mod', args.mod_bam_file, args.mod_fast5_dir)

    ivt_backbone = BackboneNetwork(args.backbone_model_path, args.extraction_layer, args.feature_width)

    unm_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)
    mod_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)

    oligo_ref_generator = OligoReferenceGenerator(ligation_ref_file=args.ref_file)
    oligo_ref_generator.collect_motif_oligos()

    for this_motif in oligo_ref_generator.motif_oligos.keys():
        unm_container.collect_motif_nucleotides(this_motif, oligo_ref_generator, enforce_ref_5mer=args.enforce_ref_5mer)
        mod_container.collect_motif_nucleotides(this_motif, oligo_ref_generator, enforce_ref_5mer=args.enforce_ref_5mer)

        if min(len(unm_container.nucleotides[this_motif]), len(mod_container.nucleotides[this_motif]))>=args.min_coverage:
            this_motif_classifier = MotifClassifier(motif=this_motif, classifier_type=args.classifier_type, scaler=args.scaler)
            this_motif_classifier.train(unm_container.nucleotides[this_motif], mod_container.nucleotides[this_motif])
            this_motif_classifier.save(
                out_model_path=os.path.join(args.classifier_model_dir, '{}.pkl'.format(this_motif)),
                draw_prc=True
            )

if __name__ == "__main__":
    main(parser.args)