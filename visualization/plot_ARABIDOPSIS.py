import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

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

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/ARABIDOPSIS/'

sel_chr = '1'

col0_file = os.path.join(source_data_dir, 'col0/chr1/mAFiA.sites.bed')
vir1_file = os.path.join(source_data_dir, 'vir1/chr1/mAFiA.sites.bed')
miclip_file = os.path.join(source_data_dir, 'parker_miclip_sites.tds')

img_out = '/home/adrian/NCOMMS_revision/images/ARABIDOPSIS'
os.makedirs(img_out, exist_ok=True)

df_col0 = pd.read_csv(col0_file, sep='\t', dtype={'chrom': str})
df_vir1 = pd.read_csv(vir1_file, sep='\t', dtype={'chrom': str})
df_merged = pd.merge(df_col0, df_vir1, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_col0', '_vir1'])

num_bins = 20
vmax = 6
ticks = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)
counts, bin_x, bin_y = np.histogram2d(
    df_merged['modRatio_vir1'], df_merged['modRatio_col0'],
    bins=[num_bins, num_bins], range=[[0, 100], [0, 100]],
)
counts_log1p = np.log10(counts + 1)

# fig_scatter = plt.figure(figsize=(3*cm, 3*cm))
# plt.scatter(df_merged['modRatio_col0'], df_merged['modRatio_vir1'], s=0.5, alpha=0.5)

fig = plt.figure(figsize=(4*cm, 4.5*cm))
ax = fig.add_subplot()
im = ax.imshow(counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
ax.set_xticks(np.linspace(0, num_bins, 5)-0.5, ticks)
ax.set_yticks(np.linspace(0, num_bins, 5)-0.5, ticks)
cbar = fig.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal', location='top')
cbar.set_ticks(np.linspace(0, vmax, 3))
ax.set_xlabel('$S_{col0}$')
ax.set_ylabel('$S_{vir1}$')
fig.savefig(os.path.join(img_out, f'hist2d_col0_vir1.{FMT}'), **fig_kwargs)

########################################################################################################################
### compare to miclip ##################################################################################################
########################################################################################################################
# from matplotlib_venn import venn2
# import pysam
# from Bio import SeqIO
# from tqdm import tqdm
#
# ref_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/reference/TAIR10_chr_all.fasta'
# ref = {}
# for record in SeqIO.parse(ref_file, 'fasta'):
#     ref[record.id] = record.seq
#

# coverage = []
# bam_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Arabidopsis_thaliana/col0/chr1/mAFiA.reads.bam'
# with pysam.AlignmentFile(bam_file, 'rb') as bam:
#     for _, row in tqdm(df_miclip_chrom.iterrows()):
#         this_cov = -1
#         for this_col in bam.pileup(row['chrom'], row['chromStart'], row['chromEnd'], truncate=True):
#             if this_col.reference_pos==row['chromStart']:
#                 this_cov = this_col.n
#         coverage.append(this_cov)
# df_miclip_chrom['coverage'] = coverage
# df_miclip_chrom.sort_values(by='chromStart', ignore_index=True, inplace=True)

########################################################################################################################
### plot miclip peaks ##################################################################################################
########################################################################################################################
miclip_start = 3122256
miclip_end = 3122304

df_miclip = pd.read_csv(miclip_file, sep='\t',
                        usecols=[0, 1, 2, 5, 6],
                        names=['chrom', 'chromStart', 'chromEnd', 'strand', 'score'],
                        dtype={'chrom': str})
df_miclip_chrom = df_miclip[df_miclip['chrom']==sel_chr]
sub_df_miclip = df_miclip_chrom[(df_miclip_chrom['chromStart']>=miclip_start) * (df_miclip_chrom['chromEnd']<=miclip_end)]

vec_x = np.arange(miclip_start, miclip_end+1)
vec_y = np.zeros_like(vec_x, dtype=np.float32)
vec_y[sub_df_miclip['chromStart'].values - miclip_start] = sub_df_miclip['score'].values

plt.figure(figsize=(6*cm, 1*cm))
plt.bar(vec_x, vec_y)
plt.xlim([miclip_start-1, miclip_end+1])
plt.ylim([0, 2])
plt.savefig(os.path.join(img_out, f'miclip_peaks_{sel_chr}_{miclip_start}_{miclip_end}.{FMT}'), **fig_kwargs)

########################################################################################################################
### count overlap with miclip ##########################################################################################
########################################################################################################################
df_col0_mod = df_col0[df_col0['modRatio']>=50]
df_col0_mod.reset_index(inplace=True, drop=True)

total_pred = len(df_col0_mod)
overlap_miclip = 0

miclip_sites = []
for _, row in df_col0_mod.iterrows():
    if np.abs(row['chromStart'] - df_miclip_chrom['chromStart'].values).min() <= 5:
        overlap_miclip += 1
        row_miclip = df_miclip_chrom.iloc[np.argmin(np.abs(row['chromStart'] - df_miclip_chrom['chromStart'].values))]
        miclip_sites.append((row_miclip['chromStart'], row_miclip['chromEnd'], row_miclip['score']))
    else:
        miclip_sites.append((-1, -1, -1))
df_miclip_sites = pd.DataFrame(miclip_sites, columns=['chromStart_miCLIP', 'chromEnd_miCLIP', 'score_miCLIP'])

df_col0_miclip_overlap = pd.concat([df_col0_mod, df_miclip_sites], axis=1)
df_col0_miclip_overlap.to_csv(os.path.join(img_out, 'col0_miclip_overlap.tsv'), sep='\t', index=False)

with open(os.path.join(img_out, 'site_count_overlap_with_miclip.tsv'), 'w') as f_out:
    f_out.write('total_mod\tmiclip_overlap\tpercentage\n')
    f_out.write(f'{total_pred}\t{overlap_miclip}\t{(overlap_miclip/total_pred):.2f}\n')

# df_mafia_miclip = pd.merge(df_col0, df_miclip_sel, on=['chrom', 'chromStart', 'chromEnd'], suffixes=['_mafia', '_miclip'])

# ref5mer = []
# for _, row in tqdm(df_miclip_chrom.iterrows()):
#     this_seq = ref[row['chrom']][(row['chromStart']-5):(row['chromStart']+6)]
#     if row['strand']=='-':
#         this_seq = this_seq.reverse_complement()
#     # np.where(np.array(this_seq) == 'A')[0]
#     ref5mer.append(str(this_seq))
# df_miclip_chrom['ref5mer'] = ref5mer

# set_col0 = set(df_col0[df_col0['modRatio']>=50]['chromStart'])
# set_miclip = set(df_miclip_chrom[df_miclip_chrom['coverage']>=50]['chromStart'])
# venn2([set_col0, set_miclip])

# plt.scatter(df_mafia_miclip['score_miclip'], df_mafia_miclip['modRatio'])
