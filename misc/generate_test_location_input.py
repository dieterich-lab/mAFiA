import pandas as pd
import os

oligo_annotation_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/oligo_reference/oligo_mod_location.bed'
ref_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/References'
test_reference_file = os.path.join(ref_dir, 'test_reference.fasta')

df_annotation = pd.read_csv(oligo_annotation_file, sep='\t')
testable_5mers = df_annotation['ref5mer'].unique()

with open(test_reference_file, 'r') as h_ref:
    ref_lines = [l.rstrip('\n') for l in h_ref.readlines()]
ref_name = ref_lines[0].lstrip('>')
ref_seq = ref_lines[1]

candidate_base = 'G'
candidate_locations = [i for i, x in enumerate(ref_seq) if x==candidate_base]

sel_loc_5mer = []
for this_loc in candidate_locations:
    if (this_loc<2) or ((len(ref_seq)-this_loc)<3):
        continue

    this_5mer = ref_seq[(this_loc-2):(this_loc+3)]
    if this_5mer in testable_5mers:
        sel_loc_5mer.append((this_loc, this_5mer))

out_tuples = [
    (
        ref_name, this_loc, this_loc+1, 'Gm', 0, '+'
    ) for this_loc, this_5mer in sel_loc_5mer
]

df_out = pd.DataFrame(out_tuples, columns=[
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
])
df_out.to_csv(os.path.join(ref_dir, 'Gm_test_sites.bed'), sep='\t', index=False)
