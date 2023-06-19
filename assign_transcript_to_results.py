import pandas as pd
import numpy as np
from Bio import SeqIO
import pysam
from tqdm import tqdm
import os

test_dataset = 'P2_WT'

results_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{test_dataset}.tsv.merged'
transcriptome_ref_file = '/home/adrian/Data/transcriptomes/GRCh38.cdna.all.fa'
transcriptome_bam_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/{test_dataset}/transcriptome_mapped.bam'
out_tsv_path = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{test_dataset}.tsv.merged.annotated'

# results_file = f'/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/train_ISA-WUE_test_{test_dataset}/res_train_ISA-WUE_test_{test_dataset}.tsv.merged'
# transcriptome_ref_file = '/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.cdna.all.fa'
# transcriptome_bam_file = f'/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/{test_dataset}/transcriptome_mapped.bam'
# out_tsv_path = f'/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/train_ISA-WUE_test_{test_dataset}/res_train_ISA-WUE_test_{test_dataset}.tsv.merged.annotated'

df_res = pd.read_csv(results_file, sep='\t')

ref = {}
for record in SeqIO.parse(transcriptome_ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

bam = pysam.AlignmentFile(transcriptome_bam_file, 'rb')
index = pysam.IndexedReads(bam)
index.build()

empty_out = 10000

df_out = pd.DataFrame()
for _, row in tqdm(df_res.iterrows()):
    read_id = row['read_id']
    read_pos = row['read_pos']
    ref_motif = row['ref_motif']
    pred_motif = row['pred_motif']

    row['annotation'] = 'unknown'

    for this_iter in index.find(read_id):
        if this_iter.flag!=0:
            continue
        ref_seq = ref[this_iter.reference_name]
        this_read_pos = read_pos
        ref_shift = 0
        this_cigar = this_iter.cigar.copy()
        while this_read_pos>0:
            next_cigar_tuple = this_cigar.pop(0)
            remaining = min(next_cigar_tuple[1], this_read_pos)
            if next_cigar_tuple[0]==0:   # match
                this_read_pos -= remaining
                ref_shift += remaining
            elif next_cigar_tuple[0]==1:   # insertion
                this_read_pos -= remaining
            elif next_cigar_tuple[0]==2:   # deletion
                ref_shift += next_cigar_tuple[1]
            elif next_cigar_tuple[0]==4:   # soft-clip
                this_read_pos -= remaining
        this_ref_pos = this_iter.reference_start + ref_shift

        # print(this_iter.flag, ref_motif, ref_seq[this_ref_pos-3:this_ref_pos+4])

        # if (this_ref_pos<this_iter.reference_end) and (ref_seq[this_ref_pos-2:this_ref_pos+3]==ref_motif):
        if (this_ref_pos<this_iter.reference_end) and (ref_motif in ref_seq[this_ref_pos - 3:this_ref_pos + 4]):
            row['annotation'] = this_iter.reference_name
            break
        # else:
        #     print('Annotation not found!')
        #     print(ref_motif, ref_seq[this_ref_pos-2:this_ref_pos+3])
    df_out = pd.concat([df_out, pd.DataFrame(row).T])

    if len(df_out)==empty_out:
        if os.path.exists(out_tsv_path):
            df_out.to_csv(out_tsv_path, sep='\t', index=False, mode='a')
        else:
            df_out.to_csv(out_tsv_path, sep='\t', index=False, mode='w')
        df_out = pd.DataFrame()

num_annotated = (df_out['annotation']!='unknown').sum()
percent_annotated = num_annotated / len(df_out) * 100
print(f'{num_annotated}/{len(df_out)} ({percent_annotated:.1f}%) annotated')