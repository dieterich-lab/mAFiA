import os
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
# fig_kwargs = dict(format=FMT, dpi=1200)
#######################################################################

def load_genome_reference(ref_file, chrs):
    print(f'Parsing genome reference {os.path.basename(ref_file)}...')
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        if record.id in chrs:
            ref[record.id] = str(record.seq)
    return ref

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/HEK293'

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
WT_bed_file = os.path.join(source_data_dir, '100_WT_0_IVT/chr6/mAFiA.sites.bed')
KO_bed_file = os.path.join(source_data_dir, 'Mettl3-KO/chr6/mAFiA.sites.bed')
img_out = '/home/adrian/NCOMMS_revision/images/METTL3'

os.makedirs(img_out, exist_ok=True)

sel_chrom = '6'
thresh_cov = 10
ref = load_genome_reference(ref_file, [sel_chrom])

########################################################################################################################
### scatter plot #######################################################################################################
########################################################################################################################
WT_bed = pd.read_csv(WT_bed_file, sep='\t', dtype={'chrom': str})
WT_bed_chrom = WT_bed[
    (WT_bed['chrom'] == sel_chrom)
    * (WT_bed['coverage'] >= thresh_cov)
    # * (df_bed['chromStart']>=sel_chromStart)
    # * (df_bed['chromEnd']<sel_chromEnd)
    ]

KO_bed = pd.read_csv(KO_bed_file, sep='\t', dtype={'chrom': str})
KO_bed_chrom = KO_bed[
    (KO_bed['chrom'] == sel_chrom)
    * (KO_bed['coverage'] >= thresh_cov)
    # * (df_bed['chromStart']>=sel_chromStart)
    # * (df_bed['chromEnd']<sel_chromEnd)
    ]

whitelist = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']
# blacklist = ['AGACC', 'AAACA', 'TAACA', 'TAACT', 'GAACA', 'TGACC']
blacklist = []

merged_bed = pd.merge(WT_bed_chrom, KO_bed_chrom, how='inner', on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'], suffixes=['_WT', '_KO'])
merged_bed_sel = merged_bed[
    ~merged_bed['ref5mer'].isin(blacklist)
    * merged_bed['ref5mer'].isin(whitelist)
]

num_bins = 20
ticks = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)

# vmax = 30
vmax = 2

fig = plt.figure(figsize=(3.5*cm, 4*cm))
# plt.scatter(merged_bed['modRatio_WT'], merged_bed['modRatio_KO'], s=0.2, alpha=0.5)
# plt.scatter(merged_bed_sel['modRatio_WT'], merged_bed_sel['modRatio_KO'], s=0.2, alpha=0.5)
# plt.plot(np.arange(0, 100), np.arange(0, 100), c='r', linewidth=0.2, alpha=0.5)
counts, bin_x, bin_y = np.histogram2d(
    merged_bed_sel['modRatio_KO'], merged_bed_sel['modRatio_WT'],
    bins=[num_bins, num_bins], range=[[0, 100], [0, 100]],
)
counts_log1p = np.log10(1+counts)
ax = fig.add_subplot()
# im = ax.imshow(counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
im = ax.imshow(counts_log1p, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
ax.axhline(y=np.where(bin_y==50)[0][0]-0.5, c='r', linestyle='--')
ax.set_xticks(np.linspace(0, num_bins, 5)-0.5, ticks)
ax.set_yticks(np.linspace(0, num_bins, 5)-0.5, ticks)
cbar = fig.colorbar(im, orientation='horizontal', location='top', fraction=0.046, pad=0.04)
cbar.set_ticks(np.linspace(0, vmax, 3))
fig.savefig(os.path.join(img_out, f'hist2d_WT_vs_KO.{FMT}'), **fig_kwargs)


########################################################################################################################
### whole-chromosome profile ###########################################################################################
########################################################################################################################
N_BINS = 10000
ref_len = len(ref[sel_chrom])

def calc_profile(in_pos, in_mod_ratios, start, end, num_bins=N_BINS):
    break_pts = np.linspace(start, end, num_bins+1)
    site_binIndex_modRatio = []
    for this_start, this_mod_ratio in zip(in_pos, in_mod_ratios):
        this_bin_ind = np.where(this_start>=break_pts)[0][-1]
        site_binIndex_modRatio.append((this_bin_ind, this_mod_ratio))
    return site_binIndex_modRatio

def calc_avg_profile(in_bin_stoichios, num_bins=N_BINS):
    binned_avg_stoichio = np.zeros(num_bins)
    for i in range(num_bins):
        all_stoichios = [stoichio for (bin, stoichio) in in_bin_stoichios if bin == i]
        if len(all_stoichios):
            # binned_avg_stoichio[i] = np.max(all_stoichios)
            # binned_avg_stoichio[i] = np.mean(all_stoichios)
            binned_avg_stoichio[i] = np.median(all_stoichios)
        else:
            binned_avg_stoichio[i] = 0
    return binned_avg_stoichio

def plot_chromosome_profile(pos, mod_ratios, plot_name, xticks=False):
    avg_profile = calc_avg_profile(calc_profile(pos, mod_ratios, 0, ref_len))

    plt.figure(figsize=(3*cm, 1.5*cm))
    plt.plot(avg_profile, linewidth=0.2)
    plt.axhline(y=50, c='r', linestyle='--', linewidth=0.2, alpha=0.5)
    plt.xlim([0, N_BINS])
    plt.ylim([0, 100])
    if xticks:
        plt.xticks(np.linspace(0.5, N_BINS-0.5, 4), [np.format_float_scientific(x, precision=1) for x in np.linspace(1, ref_len, 4)], rotation=90)
        # plt.xlabel('Chromosome Position')
    else:
        plt.xticks(np.linspace(0.5, N_BINS - 0.5, 4), [])
    # plt.ylabel('Median S')
    plt.savefig(os.path.join(img_out, plot_name), **fig_kwargs)

plot_chromosome_profile(merged_bed_sel['chromStart'].values, merged_bed_sel['modRatio_WT'].values, plot_name=f'avg_profile_WT_chr{sel_chrom}.{FMT}', xticks=False)
plot_chromosome_profile(merged_bed_sel['chromStart'].values, merged_bed_sel['modRatio_KO'].values, plot_name=f'avg_profile_METTL3_KO_chr{sel_chrom}.{FMT}', xticks=True)

########################################################################################################################
### S profile along transcript #########################################################################################
########################################################################################################################
# def parse_stoichiometry_along_gene(ranges, df_stoichio):
#     site_bin_stoichio = {}
#     for gene_name, gene_start, gene_end, gene_strand in ranges:
#         sub_df_stoichio = df_stoichio[
#             (df_stoichio['chromStart']>=gene_start) *
#             (df_stoichio['chromEnd']<gene_end) *
#             (df_stoichio['strand']==gene_strand)
#         ]
#         if len(sub_df_stoichio):
#             site_bin_stoichio[gene_name] = (calc_profile(sub_df_stoichio, gene_start, gene_end, gene_strand))
#     return site_bin_stoichio
#
#
#
# annot_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.bed'
# df_annot = pd.read_csv(annot_file, sep='\t',
#                        usecols=list(range(6)) + [12],
#                        names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'others'],
#                        dtype={'chrom': str})
#
#
# geneIDs = []
# geneNames = []
# for this_info in df_annot['others'].values:
#     splitted = this_info.split(';')
#     geneIDs.extend([this_split.split('=')[1] for this_split in splitted if this_split.split('=')[0] == 'geneID'])
#     geneNames.extend([this_split.split('=')[1] for this_split in splitted if this_split.split('=')[0] == 'gene_name'])
#
# df_annot['geneID'] = geneIDs
# df_annot['geneName'] = geneNames
# df_annot_sel = df_annot[df_annot['chrom']==sel_chrom]
#
# gene_ranges = []
# for this_gene in df_annot_sel['geneName'].unique():
#     sub_df = df_annot_sel[df_annot_sel['geneName']==this_gene]
#     gene_ranges.append((this_gene, sub_df['chromStart'].min(), sub_df['chromEnd'].max(), sub_df['strand'].unique()[0]))
#
# gene_bin_stoichios = parse_stoichiometry_along_gene(gene_ranges, df_bed_sel)
# avg_profile = calc_avg_profile(bin_stoichios)

########################################################################################################################
### parse transcript mapping ###########################################################################################
########################################################################################################################
# transcript_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/100_WT_0_IVT/transcriptome_filtered_q50.bam'
# transcriptome_file = '/home/adrian/Data/GRCh38_102/GRCh38.cdna.all.fa'
# transcriptome = {}
# for record in SeqIO.parse(transcriptome_file, 'fasta'):
#     transcriptome[record.id] = str(record.seq)
#
# def get_transcript_bin_stoichio(in_read, readPos_modProbs):
#     dict_aligned_pos = {read_pos: ref_pos for read_pos, ref_pos in in_read.get_aligned_pairs() if ref_pos is not None}
#     ref_len = len(transcriptome[in_read.reference_name])
#     bin_modProbs = []
#     for pos, prob in readPos_modProbs:
#         if pos in dict_aligned_pos.keys():
#             trans_pos = dict_aligned_pos[pos]
#             bin_ind = int((trans_pos / ref_len) * 100)
#             bin_modProbs.append((bin_ind, prob/255.0))
#     return bin_modProbs
#
#
# bam_file = '/home/adrian/NCOMMS_revision/source_data/METTL3/100WT/merged.mAFiA.reads.bam'
#
# transcript_bin_stoichio = {}
# with pysam.AlignmentFile(bam_file, 'rb') as bam_genome:
#     with pysam.AlignmentFile(transcript_file, 'rb') as bam_transcriptome:
#         for read in tqdm(bam_genome.fetch(sel_chrom)):
#             dict_mod = read.modified_bases_forward
#             if dict_mod is not None:
#                 mod_bases = dict_mod.get(('N', 0, 21891), [])
#                 if len(mod_bases):
#                     # print(mod_bases)
#                     for read_trans in bam_transcriptome.fetch():
#                         if read_trans.query_name==read.query_name:
#                             # print(read_trans, mod_bases)
#                             if read_trans.reference_name not in transcript_bin_stoichio.keys():
#                                 transcript_bin_stoichio[read_trans.reference_name] = []
#                             transcript_bin_stoichio[read_trans.reference_name].extend(get_transcript_bin_stoichio(read_trans, mod_bases))
#
#
# avg_profile = calc_avg_profile([vv for k, v in transcript_bin_stoichio.items() for vv in v])