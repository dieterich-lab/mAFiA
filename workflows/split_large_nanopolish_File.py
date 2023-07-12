import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--outfile_prefix')
parser.add_argument('--num_reads_per_file', type=int, default=400000)
args = parser.parse_args()
infile = args.infile
outfile_prefix = args.outfile_prefix
num_reads_per_file = args.num_reads_per_file

# infile = '/scratch/achan/nanopolish/100_WT_0_IVT/eventalign.txt'
# outfile_prefix = '/scratch/achan/CHEUI/100_WT_0_IVT/eventalign_part'
# num_reads_per_file = 400000

df_iterator = pd.read_csv(
    infile,
    sep='\t',
    chunksize=1E6
)

total_num_reads_written = 0
num_reads_written = 0
outfile_part = total_num_reads_written // num_reads_per_file
outfile = f'{outfile_prefix}{outfile_part:02}.txt'
df_remains = pd.DataFrame()

for df_chunk in df_iterator:
    df_continued = pd.concat([df_remains, df_chunk])
    unique_reads = df_continued['read_index'].unique()
    complete_reads = unique_reads[:-1]

    df_out = df_continued[df_continued['read_index'].isin(complete_reads)]
    df_remains = df_continued[~df_continued['read_index'].isin(complete_reads)]

    if num_reads_written==0:
        mode = 'w'
        header = True
    else:
        mode = 'a'
        header = False

    df_out.to_csv(
        outfile,
        sep='\t',
        index=False,
        header=header,
        mode=mode
    )
    num_reads_written += len(complete_reads)

    if num_reads_written>=num_reads_per_file:
        total_num_reads_written += num_reads_written
        print(f'Total num. reads written: {total_num_reads_written}')
        num_reads_written = 0
        outfile_part = total_num_reads_written // num_reads_per_file
        print(f'Staring part{outfile_part}')
        outfile = f'{outfile_prefix}{outfile_part:02}.txt'

df_remains.to_csv(
    outfile,
    sep='\t',
    index=False,
    header=False,
    mode='a'
)
num_reads_written += len(df_remains['read_index'].unique())
total_num_reads_written += num_reads_written

print(f'Total num. reads written: {total_num_reads_written}')
