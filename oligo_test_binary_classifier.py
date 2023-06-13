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

def task(in_queue, ind, containers, backbones, in_args):
    print(f"Task {containers[ind].name} started ...", flush=True)
    containers[ind].collect_features_from_reads(backbones[ind], in_args.max_num_reads)
    in_queue.put(containers[ind])
    print(f"Task {containers[ind].name} terminated with {len(containers[ind].read_bases_features)} features", flush=True)

if __name__ == "__main__":
    test_container = Oligo_Data_Container('test', args.test_bam_file, args.test_fast5_dir)
    test_container.build_dict_read_ref()

    queue = Queue()

    ivt_backbones = [Backbone_Network(args.backbone_model_path, args.extraction_layer, args.feature_width) for i in range(args.num_processes)]
    daughter_containers = test_container.get_split_containers(args.num_processes)
    processes = [Process(target=task, args=(queue, i, daughter_containers, ivt_backbones, args,)) for i in range(args.num_processes)]

    for process in processes:
        process.start()
    women_containers = []
    while not queue.empty():
        print(queue)
        women_containers.append(queue.get())
    for process in processes:
        process.join()
    print('Processes finished')
    print('{} daughters:'.format(len(daughter_containers)))
    for container in daughter_containers:
        print('{}: {}'.format(container.name, len(container.read_bases_features)))
    print('{} women:'.format(len(women_containers)))
    for container in women_containers:
        print('{}: {}'.format(container.name, len(container.read_bases_features)))
    test_container.merge_basecalls_features(daughter_containers)

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