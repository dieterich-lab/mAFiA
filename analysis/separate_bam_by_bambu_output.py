import pandas as pd
import pysam
import os
from collections import Counter

ds = 'SHAM_day56'
bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/{ds}/Fhl1.mAFiA.reads.bam'
bambu_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/{ds}/bambu'
tx_map = os.path.join(bambu_dir, 'Fhl1.tx_map.tsv')
tx_id = os.path.join(bambu_dir, 'Fhl1.tx_id.tsv')
df_tx_map = pd.read_csv(tx_map, sep='\t').convert_dtypes()
# df_tx_map = pd.read_csv(tx_map, sep='\t')
tx_ids = pd.read_csv(tx_id, sep='\t', header=None).values[:, 0]

# read_txId = {this_txId: [] for this_txId in tx_ids}
# read_txId['mixed'] = []
# read_txId['unclassified'] = []

readId_txId = {}
for _, this_row in df_tx_map.iterrows():
    this_readId = this_row['readId']
    if pd.isna(this_row['equalMatches']):
        if pd.isna(this_row['compatibleMatches']):
            # read_txId['unclassified'].append(this_readId)
            readId_txId[this_readId] = 'unclassified'
        elif str(this_row['compatibleMatches']).isnumeric():
            this_read_txId = tx_ids[int(this_row['compatibleMatches'])-1]
            # read_txId[this_read_txId].append(this_readId)
            readId_txId[this_readId] = this_read_txId
        elif str(this_row['compatibleMatches'])[0] == 'c':
            # read_txId['mixed'].append(this_readId)
            readId_txId[this_readId] = 'mixed'
    elif str(this_row['equalMatches']).isnumeric():
        this_read_txId = tx_ids[int(this_row['equalMatches'])-1]
        # read_txId[this_read_txId].append(this_readId)
        readId_txId[this_readId] = this_read_txId
    elif str(this_row['equalMatches'])[0] == 'c':
        # read_txId['mixed'].append(this_readId)
        readId_txId[this_readId] = 'mixed'

tx_counts = Counter(readId_txId.values()).most_common()
df_tx_counts = pd.DataFrame(tx_counts, columns=['transcript', 'counts'])
df_tx_counts.to_csv(os.path.join(bambu_dir, 'Fhl1_tx_counts.tsv'), sep='\t')
with pysam.AlignmentFile(bam_file, 'rb') as bam_in:
    txId_reads = {tx: [] for tx, counts in tx_counts}
    for this_read in bam_in.fetch():
        if this_read.query_name in readId_txId.keys():
            txId_reads[readId_txId[this_read.query_name]].append(this_read)
    for txId, reads in txId_reads.items():
        bam_out_file = os.path.join(bambu_dir, f'Fhl1.{txId}.bam')
        with pysam.AlignmentFile(bam_out_file, 'wb', template=bam_in) as bam_out:
            for this_read in reads:
                bam_out.write(this_read)
        pysam.index(bam_out_file)

