import os
import random

# SEED = 0
# random.seed(SEED)

ref_file = '/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta'
with open(ref_file, 'r') as f:
    all_lines = f.readlines()
ref_id_seq = {}
i = 0
while i < len(all_lines):
    if all_lines[i][0]=='>':
        id = all_lines[i].lstrip('>').rstrip('\n')
        if all_lines[i+1][0]!='>':
            seq = all_lines[i+1].rstrip('\n')
            ref_id_seq[id] = seq
            i += 2

l_margin = 23
r_margin = 33

for seed in range(100):
    random.seed(seed)

    rand_ref_id_seq = {}
    for id, seq in ref_id_seq.items():
        list_mid_seq = list(seq[l_margin:-r_margin])
        random.shuffle(list_mid_seq)
        rand_mid_seq = ''.join(list_mid_seq)
        rand_ref_id_seq[id] = seq[:l_margin] + rand_mid_seq + seq[-r_margin:]

    rand_ref_file = '/home/adrian/Data/tRNA_Berlin/newBatchDec2022_Spombe/rand_ref/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta' + '.rand' + str(seed)
    if not os.path.exists(os.path.dirname(rand_ref_file)):
        os.makedirs(os.path.dirname(rand_ref_file), exist_ok=True)
    with open(rand_ref_file, 'w') as f:
        for id, rand_seq in rand_ref_id_seq.items():
            f.write('>' + id + '\n')
            f.write(rand_seq + '\n')