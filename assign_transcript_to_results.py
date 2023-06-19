import pandas as pd
import numpy as np
from Bio import SeqIO
import pysam
from tqdm import tqdm

test_dataset = 'P2_WT'

# results_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{test_dataset}.tsv.merged'
# transcriptome_ref_file = '/home/adrian/Data/transcriptomes/GRCh38.cdna.all.fa'
# transcriptome_bam_file = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/{test_dataset}/transcriptome_mapped.bam'
# out_tsv_path = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{test_dataset}.tsv.merged.annotated'

results_file = f'/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{test_dataset}/res_train_ISA-WUE_test_{test_dataset}.tsv.merged'
transcriptome_ref_file = '/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.cdna.all.fa'
transcriptome_bam_file = f'/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/{test_dataset}/transcriptome_mapped.bam'
out_tsv_path = f'/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{test_dataset}.tsv.merged.annotated'

df_res = pd.read_csv(results_file, sep='\t')

ref = {}
for record in SeqIO.parse(transcriptome_ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

bam = pysam.AlignmentFile(transcriptome_bam_file, 'rb')
index = pysam.IndexedReads(bam)
index.build()

df_annotated = pd.DataFrame()
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

        if (this_ref_pos<this_iter.reference_end) and (ref_seq[this_ref_pos-2:this_ref_pos+3]==ref_motif):
            row['annotation'] = this_iter.reference_name
            break
    df_annotated = pd.concat([df_annotated, pd.DataFrame(row).T])
df_annotated.to_csv(out_tsv_path, sep='\t', index=False)

num_annotated = (df_annotated['annotation']!='unknown').sum()
percent_annotated = num_annotated / len(df_annotated) * 100
print(f'{num_annotated}/{len(df_annotated)} ({percent_annotated:.1f}%) annotated')

### analyze isoforms ###
# import os
# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt
# img_out = '/home/adrian/img_out/MAFIA/isoforms'
#
# min_coverage = 10
#
# isoform_mod_ratios = {}
# for this_index in tqdm(df_annotated['index'].unique()):
#     sub_df = df_annotated[df_annotated['index']==this_index]
#     sub_df = sub_df[sub_df['annotation']!='unknown']
#
#     isoform_read_nums = [(sub_df['annotation']==this_anno).sum() for this_anno in sub_df['annotation'].unique()]
#     if np.sum(np.array(isoform_read_nums)>=min_coverage)<2:
#         continue
#
#     isoform_mod_ratios[this_index] = []
#
#     for this_anno in sub_df['annotation'].unique():
#         df_anno = sub_df[sub_df['annotation']==this_anno]
#         if len(df_anno)<min_coverage:
#             continue
#         anno_mod_ratio = np.mean(df_anno['mod_prob']>=0.5)
#         isoform_mod_ratios[this_index].append((this_anno, anno_mod_ratio, len(df_anno)))
#
#     anno_ratios = [tup[1] for tup in isoform_mod_ratios[this_index]]
#     if (np.max(anno_ratios)-np.min(anno_ratios))>=0.5:
#         plt.figure(figsize=(5, 5))
#         for (this_anno, this_mod_ratio, this_num_read) in isoform_mod_ratios[this_index]:
#             plt.hist(sub_df[sub_df['annotation']==this_anno]['mod_prob'].values, range=[0, 1], bins=50, label=f"{this_anno}, {this_mod_ratio:.2f}, {this_num_read}")
#         plt.legend(loc='upper right')
#         plt.title(f"GLORI {sub_df['Ratio'].values[0]:.2f}")
#         plt.savefig(os.path.join(img_out, f"site{sub_df['index'].values[0]}.png"), bbox_inches='tight')
#         plt.close('all')