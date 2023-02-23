import os
from Bio import SeqIO
import pysam
import numpy as np
from collections import Counter
import re
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

prj_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian'
# dataset = 'A_RTA'
dataset = 'm6A_RTA'
ref_file = os.path.join(prj_dir, 'splint_variations_max_blocks_8.fasta')
fasta_file = os.path.join(prj_dir, '{}.fasta'.format(dataset))
bam_file = os.path.join(prj_dir, '{}_sorted.bam'.format(dataset))

img_out = '/home/adrian/img_out/A_m6A_RTA'
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

sequences = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    if len(record.seq)==0:
        print(record.id)
        continue
    sequences.append(record.seq)

seq_len = [len(seq) for seq in sequences]

ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = record.seq

### check bam file ###
bam = pysam.AlignmentFile(bam_file, 'rb')
ref_names = []
cs_strings = []
bam_seq_len = []
aligned_whole_ref_len = []
for record in bam.fetch():
    if record.flag==0:
        ref_names.append(record.reference_name)
        cs_strings.append(record.cigarstring)
        bam_seq_len.append(len(record.seq))
        aligned_whole_ref_len.append((record.reference_length, len(ref[record.reference_name])))

### parse cs strings ###
clip_len = []
all_ins_len = []
all_del_len = []
for cs in cs_strings:
    res = re.split('(\d+)', cs)[1:]
    cs_pairs = list(zip(res[0::2], res[1::2]))
    s_len = 0
    ins_len = []
    del_len = []
    for tup in cs_pairs:
        if tup[1]=='S':
            s_len += int(tup[0])
        elif tup[1]=='I':
            ins_len.append(int(tup[0]))
        elif tup[1]=='D':
            del_len.append(int(tup[0]))
    clip_len.append(s_len)
    all_ins_len.append(ins_len)
    all_del_len.append(del_len)

### max indel lens ###
max_ins_len = [np.max(sublist) if len(sublist)>0 else 0 for sublist in all_ins_len]
max_del_len = [np.max(sublist) if len(sublist)>0 else 0 for sublist in all_del_len]
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.hist(max_ins_len, bins=100, range=[0, 100], log=True)
plt.title('Max. insertion length in read', fontsize=15)
plt.xlabel('NT', fontsize=10)
plt.ylabel('Count', fontsize=10)
plt.subplot(1, 2, 2)
plt.hist(max_del_len, bins=100, range=[0, 100], log=True)
plt.title('Max. deletion length in read', fontsize=15)
plt.xlabel('NT', fontsize=10)
# plt.ylabel('Count', fontsize=10)
plt.suptitle(dataset, fontsize=20)
plt.savefig(os.path.join(img_out, 'hist_max_indel_per_read_{}.png'.format(dataset)), bbox_inches='tight')
plt.close('all')

plt.figure(figsize=(6, 6))
plt.plot(max_ins_len, max_del_len, '.')
plt.xlabel('Max. ins. len', fontsize=15)
plt.ylabel('Max. del. len', fontsize=15)
plt.xlim([0, 100])
plt.ylim([0, 100])
plt.suptitle(dataset, fontsize=20)
plt.savefig(os.path.join(img_out, 'scatter_max_indel_per_read_{}.png'.format(dataset)), bbox_inches='tight')
plt.close('all')

indel_thresh = 10
survival_rate = np.sum((np.array(max_ins_len)<10) * (np.array(max_del_len)<10)) / len(max_ins_len)
print('Survival rate after indel cut at {}: {:.2f}'.format(indel_thresh, survival_rate))

### block nums ###
# block_names = [name.split('_')[0] for name in ref_names]
# block_counts = Counter(block_names)
# block_num = []
# counts = []
# for k, v in block_counts.items():
#     block_num.append(int(k.lstrip('block')))
#     counts.append(v)

### plot histograms ###
# plt.figure(figsize=(10, 5))
# plt.subplot(1, 2, 1)
# plt.hist(seq_len, bins=50, range=[0, 500])
# plt.xlabel('Seq. length', fontsize=10)
# plt.ylabel('Counts', fontsize=10)
# plt.title('{} passed reads'.format(len(seq_len)), fontsize=15)
# plt.subplot(1, 2, 2)
# plt.bar(block_num, counts, width=0.5)
# plt.xlabel('Num. blocks', fontsize=10)
# # plt.ylabel('Counts', fontsize=10)
# plt.title('{} mapped reads'.format(len(block_names)), fontsize=15)
# plt.suptitle(dataset, fontsize=20)

# plt.savefig(os.path.join(img_out, '{}_read_stats.png'.format(dataset)), bbox_inches='tight')
# plt.close('all')