import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import os
from Bio import SeqIO
from Bio.Seq import Seq

thresh_modRatio = 50

# sample = 'HEK293T'
# bid_file = '/home/adrian/Data/BID_seq/BID_seq_HEK293T.bed'

sample = 'mouse heart'
bid_file = '/home/adrian/Data/BID_seq/BID_seq_mouse_heart.bed'

df_bid = pd.read_csv(
    bid_file,
    sep='\t',
    dtype={'chrom': str}
)
df_bid.rename(columns={'score': 'modRatio'}, inplace=True)

img_out = '/home/adrian/img_out/BID-Seq_motifs'
os.makedirs(img_out, exist_ok=True)

df_bid_highmod = df_bid[df_bid['modRatio']>=thresh_modRatio]
counts = Counter(df_bid_highmod['ref5mer']).most_common()
# total_count = sum([this_count[1] for this_count in counts])
# percentage = [(tup[0], tup[1]/total_count) for tup in counts]

ranked_5mers = [tup[0] for tup in counts]
modRatios = [df_bid_highmod[df_bid_highmod['ref5mer']==this_5mer]['modRatio'].values for this_5mer in ranked_5mers]

plt.figure(figsize=(6, 20))
sns.stripplot(df_bid_highmod, y='ref5mer', x='modRatio', order=ranked_5mers, hue='ref5mer', legend=False, palette='deep')
plt.xlabel('Mod. Ratio')
plt.ylabel('Motif')
plt.title(f'Bid-Seq, {sample}')
plt.savefig(os.path.join(img_out, f"{sample.replace(' ', '_')}.png"), bbox_inches='tight')

ref_file = '/home/adrian/Data/GRCm38_102/GRCm38_102.fa'

ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

oligo_list = '/home/adrian/Data/BID_seq/motifs_mouse_heart.fasta'

extended_5mers = ranked_5mers + ['TGTAG']

span = 10

with open(oligo_list, 'w') as f_out:
    for this_motif in extended_5mers:
        sub_df = df_bid[df_bid['ref5mer']==this_motif]
        freq50 = (sub_df['modRatio']>=thresh_modRatio).sum()
        row = sub_df.iloc[sub_df['modRatio'].argmax()]

        begin = row['chromStart'] - span
        end = row['chromStart'] + span + 1
        seq = ref[row['chrom']][begin:end]
        if row['strand']=='-':
            seq = str(Seq(seq).reverse_complement())

        if row['ref5mer']!=seq[span-2:span+3]:
            print(f"Warning: {row['ref5mer']} =/= {seq[span-2:span+3]}")

        if row['strand']=='+':
            seq_name = f"chr{row['chrom']}:{begin}-{end}"
        else:
            seq_name = f"chr{row['chrom']}:{end}-{begin}"

        f_out.write('>' + seq_name + f" | {row['ref5mer']}_S50_{freq50}" + '\n')
        f_out.write(seq + '\n')

