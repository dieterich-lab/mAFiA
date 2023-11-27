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
bam_file = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE/merged.mAFiA.reads.bam'
img_out = '/home/adrian/NCOMMS_revision/images/GENE_PROFILE'

########################################################################################################################
### histogram of p(m6A) at single sites ################################################################################
########################################################################################################################
sel_chrom = '6'
sel_chromStart = 31815543
sel_chromEnd = 31817946

ref = load_genome_reference(ref_file, [sel_chrom])

ref_pos_mod_probs = {}
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for col in bam.pileup(sel_chrom, sel_chromStart, sel_chromEnd, truncate=True):
        if ref[sel_chrom][col.pos]=='A':
            # print("\ncoverage at base %s = %s" % (col.pos, col.n))

            all_mod_probs = []
            for pileupread in col.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # print('\tbase in read %s = %s' %
                    #       (pileupread.alignment.query_name,
                    #        pileupread.alignment.query_sequence[pileupread.query_position]))

                    for mod_pos, mod_probs in pileupread.alignment.modified_bases_forward.get(('N', 0, 21891), []):
                        if mod_pos==pileupread.query_position:
                            all_mod_probs.append(mod_probs/255.0)

            if len(all_mod_probs):
                ref_pos_mod_probs[col.pos] = all_mod_probs

### plot histogram ###
thresh_mod = 0.5
# sel_pos = [31816450, 31817110, 31817566]
for ref_pos, mod_probs in ref_pos_mod_probs.items():
    if len(mod_probs)>=400:
    # if ref_pos in sel_pos:
        mod_ratio = int(np.mean(np.array(mod_probs)>=thresh_mod) * 100)
        plt.figure(figsize=(3*cm, 2*cm))
        plt.hist(mod_probs, bins=100, range=[0, 1])
        plt.xlim([-0.01, 1.01])
        plt.axvline(x=0.5, c='r', linestyle='--', alpha=0.5, linewidth=0.5)
        plt.title(f'chr{sel_chrom}: {ref_pos}\nS={mod_ratio}%')
        plt.savefig(os.path.join(img_out, f'hist_mod_probs_chr{sel_chrom}_{ref_pos}.{FMT}'), **fig_kwargs)
plt.close('all')

########################################################################################################################
### whole-chromosome profile ###########################################################################################
########################################################################################################################
N_BINS = 10000

def calc_profile(in_df, start, end, num_bins=N_BINS):
    break_pts = np.linspace(start, end, num_bins+1)
    site_binIndex_modRatio = []
    for _, row in in_df.iterrows():
        this_start = row['chromStart']
        this_mod_ratio = row['modRatio']
        this_bin_ind = np.where(this_start>=break_pts)[0][-1]
        site_binIndex_modRatio.append((this_bin_ind, this_mod_ratio))
    return site_binIndex_modRatio

def calc_avg_profile(in_bin_stoichios, plot_name, num_bins=N_BINS):
    binned_avg_stoichio = np.zeros(num_bins)
    for i in range(num_bins):
        all_stoichios = [stoichio for (bin, stoichio) in in_bin_stoichios if bin == i]
        if len(all_stoichios):
            # binned_avg_stoichio[i] = np.max(all_stoichios)
            binned_avg_stoichio[i] = np.median(all_stoichios)
        else:
            binned_avg_stoichio[i] = 0
    return binned_avg_stoichio

def plot_chromosome_profile(bed_file, plot_name, xticks=False):
    df_bed = pd.read_csv(bed_file, sep='\t', dtype={'chrom': str})
    df_bed_chrom = df_bed[
        (df_bed['chrom']==sel_chrom)
        # * (df_bed['chromStart']>=sel_chromStart)
        # * (df_bed['chromEnd']<sel_chromEnd)
    ]
    ref_len = len(ref['6'])
    avg_profile = calc_avg_profile(calc_profile(df_bed_chrom, 0, ref_len))

    plt.figure(figsize=(3*cm, 1.5*cm))
    plt.plot(avg_profile, linewidth=0.2)
    plt.xlim([0, N_BINS])
    plt.ylim([0, 100])
    if xticks:
        plt.xticks(np.linspace(0.5, N_BINS-0.5, 4), [np.format_float_scientific(x, precision=1) for x in np.linspace(1, ref_len, 4)], rotation=90)
        # plt.xlabel('Chromosome Position')
    else:
        plt.xticks(np.linspace(0.5, N_BINS - 0.5, 4), [])
    # plt.ylabel('Median S')
    plt.savefig(os.path.join(img_out, plot_name), **fig_kwargs)

WT_bed_file = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE/merged.mAFiA.sites.bed'
plot_chromosome_profile(WT_bed_file, plot_name=f'avg_profile_WT_chr{sel_chrom}.{FMT}', xticks=False)


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