from Bio import SeqIO
import pandas as pd
import os

oligo_ref_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/reference'
reference_file = os.path.join(oligo_ref_dir, 'WUE_oligo_ref_ABBA.fasta')
out_ref_file = os.path.join(oligo_ref_dir, 'WUE_ligation_ref_ABBA.fasta')
out_bed_file = os.path.join(oligo_ref_dir, 'WUE_annotation_ABBA.bed')

all_refs = {}
with open(reference_file, 'r') as h_ref:
    for this_record in SeqIO.parse(h_ref, format='fasta'):
        all_refs[this_record.name] = str(this_record.seq)

# df_annotation = pd.DataFrame()
key_motifs = {}

unm_loc = 6
mod_loc = 19
block_size = 26
num_blocks = 20

all_ref_keys = list(all_refs.keys())
sorted(all_ref_keys)

out_refs = {}
key_motifs = {}
for this_key in all_ref_keys:
    out_refs[this_key] = str(all_refs[this_key]) * num_blocks
    key_motifs[this_key] = all_refs[this_key][unm_loc-2:unm_loc+3]

with open(out_ref_file, 'w') as h_out_ref:
    for this_key in out_refs.keys():
        h_out_ref.write(f'>{this_key}\n')
        h_out_ref.write(out_refs[this_key] + '\n')

annotations = []
for id, motif in key_motifs.items():
    unm_locations = [i*block_size+unm_loc for i in range(num_blocks)]
    mod_locations = [i*block_size+mod_loc for i in range(num_blocks)]
    for loc in unm_locations:
        annotations.append((
            id,
            loc,
            loc+1,
            'm6A',
            255,
            '+',
            0,
            0,
            motif
        ))
    for loc in mod_locations:
        annotations.append((
            id,
            loc,
            loc+1,
            'm6A',
            255,
            '+',
            0,
            100,
            motif
        ))

df_annotations = pd.DataFrame(annotations, columns=[
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'coverage',
    'modRatio',
    'ref5mer'
])
df_annotations.sort_values(by=['chromStart'], inplace=True)
df_annotations.to_csv(out_bed_file, sep='\t', index=False)