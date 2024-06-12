import pysam
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
import os

bam_file = '/home/adrian/Data/RNA004/mapped.bam'
ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = str(record.seq)
bed_file = '/home/adrian/Data/RNA004/modProbs_per_read.bed'

header = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

df_out = []
with pysam.AlignmentFile(bam_file, "rb" ) as bam:
    for read in tqdm(bam.fetch()):
        chrom = read.reference_name
        dict_aligned_pairs = {int(k): int(v) for (k, v) in read.get_aligned_pairs(matches_only=True)}
        read_name = read.query_name
        mods = read.modified_bases
        for this_mod_key in mods.keys():
            (orig_base, strand_int, mod_base) = this_mod_key
            strand = '+' if strand_int == 0 else '-'
            pos_scores = mods[this_mod_key]
            for (this_pos, score) in pos_scores:
                chromStart = dict_aligned_pairs.get(this_pos, None)
                if chromStart is not None:
                    df_out.append([
                        chrom,
                        chromStart,
                        chromStart+1,
                        '{}_{}'.format(read_name, mod_base),
                        score,
                        strand
                    ])
        if len(df_out) >= 1000000:
            df_out = pd.DataFrame(df_out)
            if not os.path.isfile(bed_file):
                df_out.to_csv(bed_file, sep='\t', header=header, index=False)
            else:
                df_out.to_csv(bed_file, sep='\t', mode='a', header=False, index=False)
            df_out = []
df_out = pd.DataFrame(df_out)
df_out.to_csv(bed_file, sep='\t', mode='a', header=False, index=False)
