# requires:
# pip install biopython==1.81
# pip install calcs

import os, argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys
sys.path.append('~/git/mAFia_dev')
from oligo.oligo_processors import Oligo_Reference_Generator, Query_Container, Local_Aligner, Splitter, Chainer, Writer

parser = argparse.ArgumentParser(description='Map oligo basecalls through the Spanish omelette method')
parser.add_argument("--ref_file", type=str)
parser.add_argument("--query_file", type=str)
parser.add_argument("--recon_ref_file", type=str)
parser.add_argument("--sam_file", type=str)
parser.add_argument("--thresh_mapq", type=int, default=50)
parser.add_argument("--homopolymer", type=int, default=0)
parser.add_argument("--debug", default=False, action="store_true")
parser.add_argument("--write_md", default=False, action="store_true")
parser.add_argument("--write_cs", default=False, action="store_true")
args = parser.parse_args()

MIN_SEGMENT_LEN = 10

def output_score_histogram(args, in_alignments, in_queries):
    out_hist_path = os.path.join(os.path.dirname(args.sam_file),
                                 'hist_spomlette_mapping_scores_q{}.png'.format(args.thresh_mapq))
    all_mapping_scores = np.array([a.mapq for a in in_alignments])
    num_pass = np.sum(all_mapping_scores >= args.thresh_mapq)
    pass_rate = int(num_pass / len(in_queries) * 100)
    plt.figure(figsize=(5, 5))
    plt.hist(all_mapping_scores, range=[0, 100], bins=100)
    plt.xlabel('Chain Mapping scores', fontsize=12)
    plt.ylabel('Counts', fontsize=12)
    plt.axvline(x=args.thresh_mapq, c='g')
    plt.title('Pass rate at Q$\geq${}\n{}/{} = {}%'.format(args.thresh_mapq, num_pass, len(in_queries), pass_rate),
              fontsize=15)
    plt.savefig(out_hist_path, bbox_inches='tight')
    plt.close('all')

def main(args):
    oligo_ref_generator = Oligo_Reference_Generator(args.ref_file)
    queries = Query_Container(args.query_file)
    local_aligner = Local_Aligner()
    splitter = Splitter(in_aligner=local_aligner, min_segment_len=MIN_SEGMENT_LEN, thresh_mapq=args.thresh_mapq, homopolymer=args.homopolymer)
    chainer = Chainer()
    writer = Writer(out_sam_file=args.sam_file, out_ligation_ref_file=args.recon_ref_file, write_md=args.write_md, write_cs=args.write_cs)

    all_alignments = []
    print('Now mapping reads...', flush=True)
    for query in tqdm(queries.get_records()):
        if len(query.seq)<MIN_SEGMENT_LEN:
            continue
        segments = splitter.split_query_into_segments(query, oligo_ref_generator)
        if len(segments)==0:
            continue
        ref_recon = oligo_ref_generator.get_ligation_reference(segments)
        full_alignment = chainer.get_recon_align_by_chain(segments, ref_recon, query)
        all_alignments.append(full_alignment)

    output_score_histogram(args, all_alignments, queries)

    ### filter alignments by score ###
    filtered_alignments = [a for a in all_alignments if a.mapq>=args.thresh_mapq]
    filtered_ref_ids = [a.target.id for a in filtered_alignments]
    writer.write_sam(filtered_alignments, filtered_ref_ids, oligo_ref_generator)

if __name__ == "__main__":
    main(args)