from Bio import SeqIO

def add_block(base_list, add_list):
    new_list = []
    for base_block in base_list:
        for add_block in add_list:
            new_list.append(base_block+add_block)
    return new_list

base_fasta = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_splint_test_w_splint.fasta'

base_sequences = []
for record in SeqIO.parse(base_fasta, 'fasta'):
    base_sequences.append(record.seq)
base_sequences = {'seq{}'.format(i) : seq for i, seq in enumerate(base_sequences)}
base_list = list(base_sequences.values())

max_blocks = 3
current_list = [base_list.copy()]
for block in range(max_blocks-1):
    current_list.append(add_block(current_list[-1], base_list))
current_list = [seq for sublist in current_list for seq in sublist]