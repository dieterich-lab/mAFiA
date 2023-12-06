import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
import re
import pysam
from tqdm import tqdm
from Bio.Seq import Seq
from collections import Counter

def correct_mod_ratio(bam_file, df_orig):
    corr_mod_ratio = []
    corr_coverage = []
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for _, row in tqdm(df_orig.iterrows()):
            pred5mers = []
            for this_col in bam.pileup(row['chrom'], row['chromStart'], row['chromEnd'], truncate=True):
                if this_col.reference_pos == row['chromStart']:
                    # total_reads = 0
                    # corr_reads = 0
                    col_mod_probs = []
                    for pileupread in this_col.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # this_pred5mer = pileupread.alignment.get_forward_sequence()[pileupread.query_position-2:pileupread.query_position+3]
                            this_pred5mer = pileupread.alignment.query_sequence[
                                            pileupread.query_position - 2:pileupread.query_position + 3]
                            if pileupread.alignment.is_reverse:
                                this_pred5mer = str(Seq(this_pred5mer).reverse_complement())

                            if (
                                    len(this_pred5mer)
                                    and (pileupread.alignment.modified_bases is not None)
                                    and len(pileupread.alignment.modified_bases)
                                    and this_pred5mer == row['ref5mer']
                            ):
                                list_mod_probs = [pair[1] for pair in
                                                  list(pileupread.alignment.modified_bases.values())[0] if
                                                  pair[0] == pileupread.query_position]
                                if len(list_mod_probs):
                                    col_mod_probs.append(list_mod_probs[0])

                            if len(this_pred5mer):
                                pred5mers.append(this_pred5mer)

                    print(row['ref5mer'], Counter(pred5mers))
                    # print(col_mod_probs)
                    if len(col_mod_probs):
                        # pos = sum([x >= 204 for x in col_mod_probs])
                        # neg = sum([x < 51 for x in col_mod_probs])
                        # corr_mod_ratio.append(np.round(pos / (pos + neg) * 100))
                        # corr_coverage.append(pos + neg)

                        corr_mod_ratio.append(np.round(np.mean([x>=128 for x in col_mod_probs])*100))
                        corr_coverage.append(len(col_mod_probs))

                        print(this_col.n, len(col_mod_probs))
                    else:
                        corr_mod_ratio.append(0)
                        corr_coverage.append(0)

    return corr_mod_ratio, corr_coverage

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
# FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################
source_data_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_v1/Saccharomyces_cerevisiae'

dict_roman_numerals = {
    'chr1': 'I',
    'chr2': 'II',
    'chr3': 'III',
    'chr4': 'IV',
    'chr5': 'V',
    'chr6': 'VI',
    'chr7': 'VII',
    'chr8': 'VIII',
    'chr9': 'IX',
    'chr10': 'X',
    'chr11': 'XI',
    'chr12': 'XII',
    'chr13': 'XIII',
    'chr14': 'XIV',
    'chr15': 'XV',
    'chr16': 'XVI',
}

# all_chrs = ['I', 'II', 'III', 'IV']
all_chrs = ['III']
# this_chr = 'III'
all_df_wts = []
all_df_kos = []
for this_chr in all_chrs:
    wt_file = os.path.join(source_data_dir, f'WT/chr{this_chr}/mAFiA.sites.bed')
    ko_file = os.path.join(source_data_dir, f'IME4_KO/chr{this_chr}/mAFiA.sites.bed')

    all_df_wts.append(pd.read_csv(wt_file, sep='\t'))
    all_df_kos.append(pd.read_csv(ko_file, sep='\t'))
df_wt = pd.concat(all_df_wts)
df_ko = pd.concat(all_df_kos)

### bam ###
wt_bam_file = os.path.join(source_data_dir, f'WT/chr{this_chr}/mAFiA.reads.bam')
ko_bam_file = os.path.join(source_data_dir, f'IME4_KO/chr{this_chr}/mAFiA.reads.bam')

wt_mod_ratio_corr, wt_coverage_corr = correct_mod_ratio(wt_bam_file, df_wt)
df_wt['modRatio_corr'] = wt_mod_ratio_corr
df_wt['coverage_corr'] = wt_coverage_corr
ko_mod_ratio_corr, ko_coverage_corr = correct_mod_ratio(ko_bam_file, df_ko)
df_ko['modRatio_corr'] = ko_mod_ratio_corr
df_ko['coverage_corr'] = ko_coverage_corr

df_merged = pd.merge(df_wt, df_ko, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_wt', '_ko'])
df_merged = df_merged[(df_merged['coverage_corr_wt']>=50) * (df_merged['coverage_corr_ko']>=50)]

### 2d density plot, WT vs KO ###
num_bins = 40
vmax = 5
ticks = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)
counts, bin_x, bin_y = np.histogram2d(
    df_merged['modRatio_corr_ko'], df_merged['modRatio_corr_wt'],
    bins=[num_bins, num_bins], range=[[0, 100], [0, 100]],
)

fig = plt.figure(figsize=(4*cm, 4.5*cm))
ax = fig.add_subplot()
im = ax.imshow(counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
ax.set_xticks(np.linspace(0, num_bins, 5)-0.5, ticks)
ax.set_yticks(np.linspace(0, num_bins, 5)-0.5, ticks)
cbar = fig.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal', location='top')
cbar.set_ticks(np.linspace(0, vmax, 3))
ax.set_xlabel('$S_{WT}$')
ax.set_ylabel('$S_{KO}$')

### compare with mazteq ###
ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Saccharomyces_cerevisiae/reference/R64-1-1_96.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    ref[record.id] = record.seq

mod_sites_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Saccharomyces_cerevisiae/1-s2.0-S0092867419306762-mmc4.xlsx'
df_mod_sites = pd.read_excel(mod_sites_file, usecols=['chr', 'start', 'end', 'id', 'strand', 'sequence', 'confGroup', 'geneSymbol'], dtype={'start': int, 'end': int})
df_mod_sites['chrom'] = [dict_roman_numerals[this_chr] for this_chr in df_mod_sites['chr']]
df_mod_sites_filtered = df_mod_sites[
    # (df_mod_sites['confGroup']>=3)
    (df_mod_sites['chrom'].isin(all_chrs))
    ]
df_mod_sites_filtered.reset_index(drop=True, inplace=True)

df_mazter = []
for ind, row in df_mod_sites_filtered.iterrows():
    # print(row['strand'])
    ref_str = ref[row['chrom']]
    if row['strand'] == '-':
        ref_str = ref_str.reverse_complement()
    seq_match = re.search(row['sequence'], str(ref_str))
    if seq_match:
        chromStart = seq_match.span()[0] + 5
        if row['strand'] == '-':
            chromStart = len(ref_str) - chromStart - 1
        ref5mer = ref[row['chrom']][chromStart-2:chromStart+3]
        if row['strand'] == '-':
            ref5mer = ref5mer.reverse_complement()
        row['chromStart'] = chromStart
        row['chromEnd'] = chromStart + 1
        row['ref5mer'] = str(ref5mer)
        df_mazter.append(row)
    else:
        print(f'{ind}: Not found')
        # row['ref5mer'] = row['chromStart']-2:
df_mazter = pd.DataFrame(df_mazter)
df_mazter.reset_index(drop=True, inplace=True)

sel_motif = ['GGACA']
df_mazter_sel = df_mazter[df_mazter['ref5mer'].isin(sel_motif)]
df_wt_sel = df_wt[df_wt['ref5mer'].isin(sel_motif)]

df_mazter_mafia = pd.merge(df_mazter_sel, df_wt_sel, on=['chrom', 'chromStart', 'chromEnd'], suffixes=['_mazter', '_mafia'])