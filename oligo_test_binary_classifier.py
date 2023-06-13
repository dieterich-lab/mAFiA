import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/MAFIA'))
from arg_parsers import Test_Args_Parser, Output_Writer
from oligo_processors import Oligo_Reference_Generator
from data_containers import Oligo_Data_Container
from feature_extractors import Backbone_Network
from feature_classifiers import load_motif_classifiers
# from multiprocessing import Process
from torch.multiprocessing import Process, Queue

parser = Test_Args_Parser()
parser.parse_and_print()
args = parser.args

def task(ind, containers, backbones, in_args):
    containers[ind].collect_features_from_reads(backbones[ind], in_args.max_num_reads)
    queue.put(containers[ind])

if __name__ == "__main__":
    test_container = Oligo_Data_Container('test', args.test_bam_file, args.test_fast5_dir)
    test_container.build_dict_read_ref()

    queue = Queue()
    ivt_backbones = [Backbone_Network(args.backbone_model_path, args.extraction_layer, args.feature_width) for i in range(args.num_processes)]
    daughter_containers = test_container.get_split_containers(args.num_processes)
    processes = [Process(target=task, args=(i, daughter_containers, ivt_backbones, args,)) for i in range(args.num_processes)]
    for proc in processes:
        proc.start()
    women_containers = []
    for proc in processes:
        women_containers.append(queue.get())
    for proc in processes:
        proc.join()
    test_container.merge_basecalls_features(women_containers)
    print(f'After merging, collected {len(test_container.read_bases_features)} read features')

    oligo_ref_generator = Oligo_Reference_Generator(ligation_ref_file=args.ref_file)
    oligo_ref_generator.collect_motif_oligos()

    motif_classifiers = load_motif_classifiers(args.classifier_model_dir)

    writer = Output_Writer(out_path=args.outfile)

    for this_motif in oligo_ref_generator.motif_oligos.keys():
        if this_motif not in motif_classifiers.keys():
            print('No classifier available for {}. Skipping...'.format(this_motif))
            continue

        test_container.collect_motif_nucleotides(this_motif, oligo_ref_generator, enforce_ref_5mer=args.enforce_ref_5mer)

        if len(test_container.nucleotides[this_motif])<args.min_coverage:
            print('Insufficient coverage {}. Skipping...'.format(len(test_container.nucleotides[this_motif])))
            continue

        _ = motif_classifiers[this_motif].test(test_container.nucleotides[this_motif])

        df_nts = test_container.flush_nts_to_dataframe()
        writer.update_df_out(df_nts)
        writer.write_df()
    print('Total number of nucleotides tested {}'.format(len(writer.df_out)), flush=True)