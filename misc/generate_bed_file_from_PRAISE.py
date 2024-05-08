import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
        ref[record.id] = str(record.seq)

in_file = '/home/adrian/Data/PRAISE/41589_2023_1304_MOESM3_ESM.xlsx'
df_in = pd.read_excel(in_file, skiprows=range(2),
                      sheet_name='Supplementary Dataset 2',
                      usecols=[
                          'rep1_deletion_ratio',
                          'rep2_deletion_ratio',
                          'chr_site'])
out_file = '/home/adrian/Data/PRAISE/PRAISE_HEK293T_span1-3.bed'

bid_file = '/home/adrian/Data/BID_seq/BID_seq_HEK293T.bed'
df_bid = pd.read_csv(bid_file, sep='\t', dtype={'chrom': str})

df_match = []
# df_consensus = []
# df_uncertain = []
for _, row in df_in.iterrows():
    if (row['chr_site']=='NONE') or (len(row['chr_site'].split('_'))>2):
        continue

    chr = row['chr_site'].split('_')[0]
    site = row['chr_site'].split('_')[-1]
    chrom = chr.lstrip('chr')

    chromStart = -1
    chromEnd = -1

    if '-' in site:
        range_start, range_end = [int(x) for x in site.split('-')]
        strand = '+' if range_end > range_start else '-'
        range_min = min(range_start, range_end)
        range_max = max(range_start, range_end)
        span = range_max - range_min + 1
        if span==2:
            if strand=='+':
                chromStart = range_max - 1
            else:
                chromStart = range_min - 1
        elif span==3:
            chromStart = range_min
        else:
            continue
    else:
        span = 1
        chromStart = int(site) - 1
        if ref[chrom][chromStart]=='T':
            strand = '+'
        elif ref[chrom][chromStart]=='A':
            strand = '-'
        else:
            continue
    ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    if strand == '-':
        ref5mer = str(Seq(ref5mer).reverse_complement())
    print(strand, span, ref5mer, )

    ### ambiguous pos ###
    # if '-' in site:
    #     flag_consensus = True
    #     ranges = [int(x) for x in site.split('-')]
    #     range_start = min(ranges)
    #     range_end = max(ranges)
    #     consensus = df_bid[(df_bid['chrom'] == chrom) * (df_bid['chromStart']>=range_start) * (df_bid['chromEnd']<=range_end)]
    #     if len(consensus)==1:
    #         # print(row)
    #         # print(consensus)
    #         # print('\n')
    #         chromStart = consensus['chromStart'].values[0]
    #         chromEnd = consensus['chromEnd'].values[0]
    #
    # else:
    #     flag_consensus = False
    #     chromEnd = int(site)
    #     chromStart = chromEnd - 1
    #     # strand = row['strand']
    #
    # # ref5mer = ref[chrom][chromStart-2:chromStart+3]
    # # if strand=='-':
    # #     ref5mer = str(Seq(ref5mer).reverse_complement())
    #
    # if chromStart<0:
    #     continue

    test_chars = ref5mer[2], ref5mer[1], ref5mer[3]
    if [test_chars[i]=='T' for i in range(span)]:
        chromEnd = chromStart + 1
        avg_deletion_ration = round((row['rep1_deletion_ratio'] + row['rep2_deletion_ratio']) * 0.5 * 100.0, 1)
        df_match.append([chrom, chromStart, chromEnd, 'psi', avg_deletion_ration, strand, ref5mer])

    # if ref5mer[2]=='T':
    #     strand = '+'
    # elif ref5mer[2]=='A':
    #     ref5mer = str(Seq(ref5mer).reverse_complement())
    #     strand = '-'
    # else:
    #     continue

    # if strand=='-':
    #     ref9mer = str(Seq(ref9mer).reverse_complement())
    #
    # span = re.search(row['Motif_1'], ref9mer).span()
    # if span[0]!=2:
    #     continue
    #     # shift = span[0] - 2
    #     # if strand=='+':
    #     #     chromStart += shift
    #     #     chromEnd += shift
    #     # else:
    #     #     chromStart -= shift
    #     #     chromEnd -= shift
    # ref5mer = ref[chrom][chromStart - 2:chromStart + 3]
    # if strand == '-':
    #     ref5mer = str(Seq(ref5mer).reverse_complement())

    # print(ref5mer==row['Motif_1'])

    # if flag_consensus:
    #     df_consensus.append([chrom, chromStart, chromEnd, 'psi', avg_deletion_ration, strand, ref5mer])
    # else:
    #     df_match.append([chrom, chromStart, chromEnd, 'psi', avg_deletion_ration, strand, ref5mer])

df_match = pd.DataFrame(df_match, columns=bed_fields)
# df_consensus = pd.DataFrame(df_consensus, columns=bed_fields)
# df_out = pd.concat([df_match, df_consensus])
df_out = df_match
df_out = pd.concat([
    df_out[df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart'], key=pd.to_numeric),
    df_out[~df_out['chrom'].str.isnumeric()].sort_values(by=['chrom', 'chromStart']),
])
df_out.to_csv(out_file, sep='\t', index=False)