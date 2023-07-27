import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
from arg_parsers import Train_Args_Parser
from oligo_processors import Oligo_Reference_Generator
from data_containers import Oligo_Data_Container
from feature_extractors import Backbone_Network
from feature_classifiers import Motif_Classifier

parser = Train_Args_Parser()
parser.parse_and_print()

def main(args):
    os.makedirs(args.classifier_model_dir, exist_ok=True)

    unm_container = Oligo_Data_Container('unm', args.unm_bam_file, args.unm_fast5_dir)
    mod_container = Oligo_Data_Container('mod', args.mod_bam_file, args.mod_fast5_dir)

    ivt_backbone = Backbone_Network(args.backbone_model_path, args.extraction_layer, args.feature_width)

    unm_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)
    mod_container.collect_features_from_reads(ivt_backbone, args.max_num_reads)

    oligo_ref_generator = Oligo_Reference_Generator(ligation_ref_file=args.ref_file)
    oligo_ref_generator.collect_motif_oligos()

    for this_motif in oligo_ref_generator.motif_oligos.keys():
        unm_container.collect_motif_nucleotides(this_motif, oligo_ref_generator, enforce_ref_5mer=args.enforce_ref_5mer)
        mod_container.collect_motif_nucleotides(this_motif, oligo_ref_generator, enforce_ref_5mer=args.enforce_ref_5mer)

        if min(len(unm_container.nucleotides[this_motif]), len(mod_container.nucleotides[this_motif]))>=args.min_coverage:
            this_motif_classifier = Motif_Classifier(motif=this_motif, classifier_type=args.classifier_type, scaler=args.scaler)
            this_motif_classifier.train(unm_container.nucleotides[this_motif], mod_container.nucleotides[this_motif])
            this_motif_classifier.save(
                out_model_path=os.path.join(args.classifier_model_dir, '{}.pkl'.format(this_motif)),
                draw_prc=True
            )

if __name__ == "__main__":
    main(parser.args)