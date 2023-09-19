from Bio import SeqIO
import pandas as pd
import os

reference_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/References/train_reference.fasta'
oligo_ref_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/oligo_reference'
os.makedirs(oligo_ref_dir, exist_ok=True)

all_refs = {}
with open(reference_file, 'r') as h_ref:
    for this_record in SeqIO.parse(h_ref, format='fasta'):
        all_refs[this_record.name] = str(this_record.seq)

# df_annotation = pd.DataFrame()
key_motifs = {}

unm_loc = 5
mod_loc = 16

all_ref_keys = list(all_refs.keys())
sorted(all_ref_keys)
for group in ['A', 'B', 'C', 'D']:
    group_keys = [k for k in all_ref_keys if k.split('-')[0]==group]
    group_size = len(group_keys)
    for start in range(1, group_size+1, 4):
        stop = start + 3

        outfile = os.path.join(oligo_ref_dir, f'oligo_{group}{start}-{stop}.fasta')

        subgroup_keys = [f'{group}-{ind:02d}' for ind in range(start, stop+1, 1)]
        with open(outfile, 'w') as h_out_ref:
            for this_key in subgroup_keys:
                out_key = 'GM-' + this_key.replace('-', '')
                h_out_ref.write(f'>{out_key}\n')
                h_out_ref.write(all_refs[this_key]+'\n')

                this_motif_unm = all_refs[this_key][(unm_loc-2):(unm_loc+3)]
                this_motif_mod = all_refs[this_key][(mod_loc-2):(mod_loc+3)]
                if this_motif_unm!=this_motif_mod:
                    print(f'Error! {this_motif_unm}=/={this_motif_mod}')
                    break

                key_motifs[out_key] = this_motif_unm

annotations = []
for id, motif in key_motifs.items():
    for loc in [unm_loc, mod_loc]:
        annotations.append((
            id,
            loc,
            loc+1,
            'Gm',
            255,
            '+',
            0,
            0 if loc==unm_loc else 100,
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
df_annotations.to_csv(os.path.join(oligo_ref_dir, 'oligo_mod_location.bed'), sep='\t', index=False)