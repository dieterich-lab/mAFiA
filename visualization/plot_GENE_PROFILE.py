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

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
annot_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.bed'
transcript_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/100_WT_0_IVT/transcriptome_filtered_q50.bam'
bam_file = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE/100WT/merged.mAFiA.reads.bam'
img_out = '/home/adrian/NCOMMS_revision/images/GENE_PROFILE'

########################################################################################################################
### histogram of p(m6A) at single sites ################################################################################
########################################################################################################################
sel_chrom = '6'
sel_chromStart = 31815543
sel_chromEnd = 31817946

ref = load_genome_reference(ref_file, [sel_chrom])

# ref_pos_mod_probs = {}
# with pysam.AlignmentFile(bam_file, 'rb') as bam:
#     for col in bam.pileup(sel_chrom, sel_chromStart, sel_chromEnd, truncate=True):
#         if ref[sel_chrom][col.pos]=='A':
#             # print("\ncoverage at base %s = %s" % (col.pos, col.n))
#
#             all_mod_probs = []
#             for pileupread in col.pileups:
#                 if not pileupread.is_del and not pileupread.is_refskip:
#                     # print('\tbase in read %s = %s' %
#                     #       (pileupread.alignment.query_name,
#                     #        pileupread.alignment.query_sequence[pileupread.query_position]))
#
#                     for mod_pos, mod_probs in pileupread.alignment.modified_bases_forward.get(('N', 0, 21891), []):
#                         if mod_pos==pileupread.query_position:
#                             all_mod_probs.append(mod_probs/255.0)
#
#             if len(all_mod_probs):
#                 ref_pos_mod_probs[col.pos] = all_mod_probs

### plot histogram ###
# thresh_mod = 0.5
# sel_pos = [31816312, 31816986, 31817582]
# for ref_pos, mod_probs in ref_pos_mod_probs.items():
#     # if len(mod_probs)>=400:
#     if ref_pos in sel_pos:
#         mod_ratio = int(np.mean(np.array(mod_probs)>=thresh_mod) * 100)
#         plt.figure(figsize=(3*cm, 2*cm))
#         plt.hist(mod_probs, bins=100, range=[0, 1])
#         plt.xlim([-0.01, 1.01])
#         plt.axvline(x=0.5, c='r', linestyle='--', alpha=0.5, linewidth=0.5)
#         plt.title(f'chr{sel_chrom}: {ref_pos}\nS={mod_ratio}%')
#         plt.savefig(os.path.join(img_out, f'hist_mod_probs_chr{sel_chrom}_{ref_pos}.{FMT}'), **fig_kwargs)
# plt.close('all')

########################################################################################################################
### filter bam #########################################################################################################
########################################################################################################################
def generate_mm_ml_tags(read_mods, mod_base='N', mod_code='21891'):
    dists = [read_mods[0][0]] + list(np.diff([mod[0] for mod in read_mods]) - 1)
    mod_probs = [mod[1] for mod in read_mods]
    strand = '+'
    mm_tag = mod_base + strand + mod_code + ',' + ','.join([str(d) for d in dists]) + ';'
    ml_tag = mod_probs
    return mm_tag, ml_tag

sel_motifs = [
    'AGACT',
    'GAACT',
    'GGACA',
    'GGACC',
    'GGACT',
    'TGACT',
]

# out_bam_file = bam_file.replace('merged', f'chr{sel_chrom}_{sel_chromStart}_{sel_chromEnd}').replace('reads.bam', 'reads.6motifs.bam')
# with pysam.AlignmentFile(bam_file, 'rb') as bam_in:
#     with pysam.AlignmentFile(out_bam_file, 'wb', template=bam_in) as bam_out:
#         for read in bam_in.fetch(sel_chrom, sel_chromStart, sel_chromEnd):
#             old_read_mods = read.modified_bases_forward.get(('N', 0, 21891), [])
#             if len(old_read_mods):
#                 dict_aligned_pos = {read_pos: ref_pos for read_pos, ref_pos in read.get_aligned_pairs() if
#                                     ref_pos is not None}
#                 new_read_mods = []
#                 for mod_pos, mod_probs in old_read_mods:
#                     motif = ref[sel_chrom][dict_aligned_pos[mod_pos]-2:dict_aligned_pos[mod_pos]+3]
#                     if motif in sel_motifs:
#                         # print(motif)
#                         new_read_mods.append((mod_pos, mod_probs))
#                 if len(new_read_mods):
#                     new_mm, new_ml = generate_mm_ml_tags(new_read_mods)
#                     read.set_tag('MM', new_mm)
#                     read.set_tag('ML', new_ml)
#                     bam_out.write(read)




########################################################################################################################
### scatter plot #######################################################################################################
########################################################################################################################
sel_chrom = '6'
thresh_cov = 50

# WT_bed_file = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE/100WT/merged.mAFiA.sites.bed'
# KO_bed_file = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE/Mettl3-KO/chr6/mAFiA.sites.bed'

WT_bed_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/100_WT_0_IVT/chr6/mAFiA.sites.bed'
KO_bed_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_replaced/Mettl3-KO/chr6/mAFiA.sites.bed'

# WT_bed_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/_DRACH_var_thresh/100_WT_0_IVT/chr6/mAFiA.sites.bed'
# KO_bed_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/_DRACH_var_thresh/Mettl3-KO/chr6/mAFiA.sites.bed'


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
    # * merged_bed['ref5mer'].isin(whitelist)
]

num_bins = 20
ticks = np.int32(np.linspace(0, num_bins, 5) * 100 / num_bins)

vmax = 6

fig = plt.figure(figsize=(3.5*cm, 4*cm))
# plt.scatter(merged_bed['modRatio_WT'], merged_bed['modRatio_KO'], s=0.2, alpha=0.5)
# plt.scatter(merged_bed_sel['modRatio_WT'], merged_bed_sel['modRatio_KO'], s=0.2, alpha=0.5)
# plt.plot(np.arange(0, 100), np.arange(0, 100), c='r', linewidth=0.2, alpha=0.5)
counts, bin_x, bin_y = np.histogram2d(
    merged_bed_sel['modRatio_KO'], merged_bed_sel['modRatio_WT'],
    bins=[num_bins, num_bins], range=[[0, 100], [0, 100]],
)
ax = fig.add_subplot()
im = ax.imshow(counts, origin='lower', cmap=mpl.cm.plasma, vmin=0, vmax=vmax)
ax.set_xticks(np.linspace(0, num_bins, 5)-0.5, ticks)
ax.set_yticks(np.linspace(0, num_bins, 5)-0.5, ticks)
cbar = fig.colorbar(im, orientation='horizontal', location='top', fraction=0.046, pad=0.04)
cbar.set_ticks(np.linspace(0, vmax, 3))
fig.savefig(os.path.join(img_out, f'hist2d_WT_vs_KO.{FMT}'), **fig_kwargs)


########################################################################################################################
### whole-chromosome profile ###########################################################################################
########################################################################################################################
N_BINS = 100000
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