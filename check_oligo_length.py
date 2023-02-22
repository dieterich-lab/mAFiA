from Bio import SeqIO
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# fasta_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/A_RTA.fasta'
fasta_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A_RTA.fasta'

sequences = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    if len(record.seq)==0:
        print(record.id)
        continue
    sequences.append(record.seq)

seq_len = [len(seq) for seq in sequences]

plt.figure(figsize=(8, 8))
plt.hist(seq_len, bins=250, range=[0, 250])