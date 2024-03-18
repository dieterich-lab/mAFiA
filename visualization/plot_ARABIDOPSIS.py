import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
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

                    print(row['ref5mer'], Counter(pred5mers).most_common())
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

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/ARABIDOPSIS/'
miclip_file = os.path.join(source_data_dir, 'parker_miclip_sites.tds')

img_out = '/home/adrian/NCOMMS_revision/images/ARABIDOPSIS'
os.makedirs(img_out, exist_ok=True)

sel_chrs = ['1', '2', '3', '4', '5']

df_col0 = []
df_vir1 = []
for this_chr in sel_chrs:
    col0_file = os.path.join(source_data_dir, f'col0/chr{this_chr}/mAFiA.sites.bed')
    vir1_file = os.path.join(source_data_dir, f'vir1/chr{this_chr}/mAFiA.sites.bed')
    col0_bam_file = os.path.join(source_data_dir, f'col0/chr{this_chr}/mAFiA.reads.bam')
    vir1_bam_file = os.path.join(source_data_dir, f'vir1/chr{this_chr}/mAFiA.reads.bam')
    df_col0.append(pd.read_csv(col0_file, sep='\t', dtype={'chrom': str}))
    df_vir1.append(pd.read_csv(vir1_file, sep='\t', dtype={'chrom': str}))
df_col0 = pd.concat(df_col0, ignore_index=True)
df_vir1 = pd.concat(df_vir1, ignore_index=True)

# col0_mod_ratio_corr, col0_coverage_corr = correct_mod_ratio(col0_bam_file, df_col0)
# df_col0['modRatio_corr'] = col0_mod_ratio_corr
# df_col0['coverage_corr'] = col0_coverage_corr
#
# vir1_mod_ratio_corr, vir1_coverage_corr = correct_mod_ratio(vir1_bam_file, df_vir1)
# df_vir1['modRatio_corr'] = vir1_mod_ratio_corr
# df_vir1['coverage_corr'] = vir1_coverage_corr

df_merged = pd.merge(df_col0, df_vir1, on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_col0', '_vir1'])
df_merged = df_merged[(df_merged['coverage_col0']>=50) * (df_merged['coverage_vir1']>=50)]

sel_motifs = [
    'GGACT',
    'GGACA',
    'GAACT',
    'AGACT',
    'GGACC',
    'TGACT'
]

# plt.figure(figsize=(10, 6))
# for ind, motif in enumerate(sel_motifs):
#     sub_df = df_merged[df_merged['ref5mer']==motif]
#     plt.subplot(2, 3, ind+1)
#     plt.scatter(sub_df['modRatio_col0'], sub_df['modRatio_vir1'], c='b', s=1, alpha=0.5)
#     # plt.scatter(sub_df['modRatio_corr_col0'], sub_df['modRatio_corr_vir1'], s=1.5, c='r', alpha=0.5)
#     plt.xlim([-1, 101])
#     plt.ylim([-1, 101])
#     plt.title(motif)

df_merged_sel = df_merged[df_merged['ref5mer'].isin(sel_motifs)]
num_bins = 20
# vmax = 2
vmax = 20
ticks = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)
counts, bin_x, bin_y = np.histogram2d(
    df_merged_sel['modRatio_vir1'], df_merged_sel['modRatio_col0'],
    bins=[num_bins, num_bins], range=[[0, 100], [0, 100]],
)
counts_log1p = np.log10(counts + 1)

# fig_scatter = plt.figure(figsize=(3*cm, 3*cm))
# plt.scatter(df_merged['modRatio_col0'], df_merged['modRatio_vir1'], s=0.5, alpha=0.5)

fig_hist2d = plt.figure(figsize=(4.5*cm, 4*cm))
ax_hist2d = fig_hist2d.add_subplot()
im = ax_hist2d.imshow(counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
# im = ax_hist2d.imshow(counts_log1p, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
# ax_hist2d.axhline(y=np.where(bin_y==50)[0][0]-0.5, c='r', linestyle='--')
ax_hist2d.plot([0, num_bins-1], [0, num_bins-1], c='r', linestyle='--', alpha=0.5)
ax_hist2d.set_xticks(np.linspace(0, num_bins, 5)-0.5, ticks)
ax_hist2d.set_yticks(np.linspace(0, num_bins, 5)-0.5, ticks)
cbar = fig_hist2d.colorbar(im, fraction=0.046, pad=0.04, orientation='horizontal', location='top')
cbar.set_ticks(np.linspace(0, vmax, 3))
# cbar.set_label('$log_{10}(1+count)$', rotation=-90, labelpad=10)
ax_hist2d.set_xlabel('$S_{col0}$')
ax_hist2d.set_ylabel('$S_{vir1}$')
fig_hist2d.savefig(os.path.join(img_out, f'hist2d_col0_vir1.{FMT}'), **fig_kwargs)

with open(os.path.join(source_data_dir, 'source_data_Figure_2h.tsv'), 'w') as fout:
    fout.write('Figure 2h\n\n')
    fout.write('\t' + 'S_vir1 \ S_col0' + '\t' + '\t'.join([str(int(x)) for x in bin_x[:-1]]) + '\n')
    for row, this_bin_y in enumerate(bin_y[:-1]):
        fout.write('\t' + str(int(this_bin_y)))
        fout.write('\t' + '\t'.join([str(int(x)) for x in counts[row]]) + '\n')

### 1D Histogram ###
# df_merged_mod = df_merged[df_merged['modRatio_col0']>=50]
# # df_merged_mod = df_merged
#
# fig_hist1d = plt.figure(figsize=(4*cm, 4*cm))
# ax_hist1d = fig_hist1d.add_subplot()
# ax_hist1d.hist(df_merged_mod['modRatio_col0'], bins=num_bins, range=[0, 100], alpha=0.5, label='col0')
# ax_hist1d.hist(df_merged_mod['modRatio_vir1'], bins=num_bins, range=[0, 100], alpha=0.5, label='vir1')
# # ax_hist1d.hist(df_merged_mod['modRatio_col0'], bins=40, range=[0, 100], alpha=0.5, label='col0', log=True)
# # ax_hist1d.hist(df_merged_mod['modRatio_vir1'], bins=40, range=[0, 100], alpha=0.5, label='vir1', log=True)
# ax_hist1d.set_xlim([-1, 101])
# ax_hist1d.legend(loc='upper right', handlelength=1)
# ax_hist1d.set_xlabel('S')
# ax_hist1d.set_ylabel('Num. Sites')
# # ax_hist1d.set_title(f'{len(df_merged_mod)} sites' + ' with $S_{col0}{\geq}$50')
# fig_hist1d.savefig(os.path.join(img_out, f'hist1d_col0_vir1.{FMT}'), **fig_kwargs)

########################################################################################################################
### compare to miclip ##################################################################################################
########################################################################################################################
df_miclip = pd.read_csv(miclip_file, sep='\t',
                        usecols=[0, 1, 2, 5, 6],
                        names=['chrom', 'chromStart', 'chromEnd', 'strand', 'score'],
                        dtype={'chrom': str, 'chromStart': int, 'chromEnd': int})

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
miclip_chr = '1'
miclip_start = 3122256
miclip_end = 3122304

df_miclip_chrom = df_miclip[df_miclip['chrom']==miclip_chr]
sub_df_miclip = df_miclip_chrom[(df_miclip_chrom['chromStart']>=miclip_start) * (df_miclip_chrom['chromEnd']<=miclip_end)]

vec_x = np.arange(miclip_start, miclip_end+1)
vec_y = np.zeros_like(vec_x, dtype=np.float32)
vec_y[sub_df_miclip['chromStart'].values - miclip_start] = sub_df_miclip['score'].values

with open(os.path.join(source_data_dir, 'source_data_Figure_2g.tsv'), 'w') as fout:
    fout.write('Figure 2g\n\n')
    fout.write('\t' + 'chr1 pos' + '\t' + '\t'.join([str(int(x)) for x in vec_x]) + '\n')
    fout.write('\t' + 'miCLiP log2FC' + '\t' + '\t'.join([str(x) for x in vec_y]) + '\n')



plt.figure(figsize=(6*cm, 1*cm))
plt.bar(vec_x, vec_y)
plt.xlim([miclip_start-1, miclip_end+1])
plt.ylim([0, 2])
plt.savefig(os.path.join(img_out, f'miclip_peaks_{miclip_chr}_{miclip_start}_{miclip_end}.{FMT}'), **fig_kwargs)

########################################################################################################################
### count overlap with miclip ##########################################################################################
########################################################################################################################
df_col0_mod = df_col0[df_col0['modRatio']>=50]
df_col0_mod.reset_index(inplace=True, drop=True)

total_pred = len(df_col0_mod)
overlap_miclip = 0

miclip_sites = []
for _, row in df_col0_mod.iterrows():
    df_miclip_chrom = df_miclip[df_miclip['chrom']==row['chrom']]
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
