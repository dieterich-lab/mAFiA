import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def add_block(base_list, add_list, base_names, add_names):
    new_list = []
    new_names = []
    for (base_block, base_name) in zip(base_list, base_names):
        for (add_block, add_name) in zip(add_list, add_names):
            new_list.append(base_block+add_block)
            new_names.append(
                'block{}_{}'.format
                    (
                        int(base_name.split('_')[0].lstrip('block'))+int(add_name.split('_')[0].lstrip('block')),
                        base_name.split('_')[1] + add_name.split('_')[1]
                    )
            )
    return new_list, new_names

max_blocks = 8
ref_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/reference'
### WUE batch 1 ###
# base_fasta = os.path.join(ref_dir, 'WUE_batch1_w_splint.fasta')
# out_fasta = os.path.join(ref_dir, 'WUE_batch1_blocks{}.fasta'.format(max_blocks))
### WUE batch 2 ###
# base_fasta = os.path.join(ref_dir, 'WUE_batch2_w_splint.fasta')
# out_fasta = os.path.join(ref_dir, 'WUE_batch2_blocks{}.fasta'.format(max_blocks))
### Isabel all 6 ###
# base_fasta = os.path.join(ref_dir, 'Top-6-GLORI_from_48-oligo-list.fasta')
# out_fasta = os.path.join(ref_dir, 'top6_random_permutation_max_blocks_{}.fasta'.format(max_blocks))
### Isabel 3+3 ###
# base_fasta = os.path.join(ref_dir, 'RL_Mix1_Mix3.fasta')
# out_fasta = os.path.join(ref_dir, 'RL_Mix1_Mix3_blocks{}.fasta'.format(max_blocks))
base_fasta = os.path.join(ref_dir, 'RL_Mix2_Mix4.fasta')
out_fasta = os.path.join(ref_dir, 'RL_Mix2_Mix4_blocks{}.fasta'.format(max_blocks))

base_sequences = []
for record in SeqIO.parse(base_fasta, 'fasta'):
    base_sequences.append(record.seq)
base_sequences = {'seq{}'.format(i) : seq for i, seq in enumerate(base_sequences)}
base_seqs = list(base_sequences.values())
base_names = ['block1_{}'.format(i) for i in range(len(base_seqs))]

current_seqs = [base_seqs.copy()]
current_names = [base_names.copy()]
for block in range(max_blocks-1):
    extra_seqs, extra_names = add_block(current_seqs[-1], base_seqs, current_names[-1], base_names)
    current_seqs.append(extra_seqs)
    current_names.append(extra_names)
current_seqs = [seq for sublist in current_seqs for seq in sublist]
current_names = [name for sublist in current_names for name in sublist]

### write out ###
seqRecords = []
for seq, name in zip(current_seqs, current_names):
    seqRecords.append(
        SeqRecord(
            seq,
            id = name,
            description = '',
            name = ''
        )
    )
with open(out_fasta, 'w') as output_handle:
    SeqIO.write(seqRecords, output_handle, "fasta")