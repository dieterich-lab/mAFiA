import os
HOME = os.path.expanduser('~')
from Bio import SeqIO
from collections import Counter
import numpy as np
import random
randSeed = 2
random.seed(randSeed)

MAX_MOTIFS = 10

def get_num_mat_from_alpha(list_bcs):
    alpha_to_num = {
        'A' : 1,
        'C' : 2,
        'G' : 3,
        'T' : 4
    }

    list_nums = [[alpha_to_num[x] for x in bc] for bc in list_bcs]
    return np.vstack(list_nums)


def get_distance_mat(mat0, mat1):
    return np.int32((mat0[:, np.newaxis, :] - mat1[np.newaxis, :, ])!=0).sum(axis=-1)

### paths ###
data_dir = os.path.join(HOME, 'Data/DRACH')
in_fasta_file = os.path.join(data_dir, 'oligoSlop8.fasta')
out_fasta_file = os.path.join(data_dir, 'chosen_DRACH_oligoSlop8_{}motifs_randSeed{}.fasta'.format(MAX_MOTIFS, randSeed))

drachs = []
barcodes = []
for record in SeqIO.parse(in_fasta_file, "fasta"):
    seq = record.seq
    if (len(seq)!=21) or ('GGGG' in seq):
        continue
    drachs.append(str(seq[8:13]))
    barcodes.append(str(seq[:8] + seq[13:]))

motif_counts = Counter(drachs).most_common()
for motif, count in motif_counts:
    print(motif, count)


### first two motifs ###
print('\nn = 0')

motif_barcode = {}
motif0 = motif_counts[0][0]
motif1 = motif_counts[1][0]
bcs0 = [bc for (drh, bc) in zip(drachs, barcodes) if drh==motif0]
bcs1 = [bc for (drh, bc) in zip(drachs, barcodes) if drh==motif1]
mat_bcs0 = get_num_mat_from_alpha(bcs0)
mat_bcs1 = get_num_mat_from_alpha(bcs1)
dist_mat = get_distance_mat(mat_bcs0, mat_bcs1)

# deterministic #
# ind0, ind1 = np.unravel_index(dist_mat.argmax(), dist_mat.shape)

# random #
candidate_ind0, candidate_ind1 = np.where(dist_mat==dist_mat.max())
r_ind = random.sample(range(len(candidate_ind0)), k=1)[0]
ind0 = candidate_ind0[r_ind]
ind1 = candidate_ind1[r_ind]

motif_barcode[motif0] = bcs0[ind0]
motif_barcode[motif1] = bcs1[ind1]

print('Chosen new barcodes:')
for bc in motif_barcode.values():
    print('\t', bc)
print('with distance', dist_mat[ind0, ind1])

### n >= 2 ###
for n in range(2, MAX_MOTIFS):
    print('\nn = {}'.format(n))
    bcs_old = list(motif_barcode.values())
    print('Existing barcodes:')
    for bc in bcs_old:
        print('\t', bc)

    motif_new = motif_counts[n][0]
    bcs_new = [bc for (drh, bc) in zip(drachs, barcodes) if drh==motif_new]

    dist_mat = get_distance_mat(get_num_mat_from_alpha(bcs_old), get_num_mat_from_alpha(bcs_new))
    ind_new = dist_mat.prod(axis=0).argmax()

    motif_barcode[motif_new] = bcs_new[ind_new]

    print('Chosen new barcode:')
    print('\t', motif_barcode[motif_new], 'with distances:', dist_mat[:, ind_new])

### summary ###
with open(out_fasta_file, "w") as fout:
    for motif in motif_barcode.keys():
        count = [v for k, v in motif_counts if k==motif][0]
        read_name = '>{}_freq{}'.format(motif, count)
        barcode = motif_barcode[motif]
        full_seq = barcode[:8] + motif + barcode[8:]

        fout.write('{}\n'.format(read_name))
        fout.write('{}\n'.format(full_seq))

for motif in motif_barcode.keys():
    if ((motif[0] in ['A', 'G', 'T']) and (motif[1] in ['A', 'G']) and (motif[2]=='A') and (motif[3]=='C') and (motif[4] in ['A', 'C', 'T'])):
        print(motif, 'DRACH')
    else:
        print(motif, 'NOT')