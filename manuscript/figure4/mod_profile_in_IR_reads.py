import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from tqdm import tqdm
import pysam
import pybedtools
from collections import Counter
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

# ds = 'TAC'
# cond = 'SHAM_merged'

# ds = 'HFpEF'
# cond = 'HFpEF_merged'

ds = 'Diet'
cond = 'WT_WD_merged'

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}
mod_colors = {
    'm6A': 'red',
    'psi': 'blue'
}
mod_labels = [(mpatches.Patch(color=this_color), this_label)
              for this_label, this_color in mod_colors.items()]

def calc_mod_level_per_read(in_reads, in_num_bases=5):
    out_read_avg_mod_level = {}
    for this_in_read in in_reads:
        out_read_avg_mod_level[this_in_read.query_name] = {}
        for mod in dict_mod_code.keys():
            mod_bases = this_in_read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])
            if len(mod_bases) >= in_num_bases:
                out_read_avg_mod_level[this_in_read.query_name][mod] = np.mean(
                    np.partition([tup[1] for tup in mod_bases], -in_num_bases)[-in_num_bases:]) / 255.0
            else:
                out_read_avg_mod_level[this_in_read.query_name][mod] = np.mean([tup[1] for tup in mod_bases]) / 255.0
    return out_read_avg_mod_level


roi_max_bin = 150
roi_bin_width = 50
roi_bin_edges = np.arange(-roi_max_bin, roi_max_bin+1, roi_bin_width)
x_zero_ind = np.where(np.array(roi_bin_edges) == 0)[0][0]
x_zero = x_zero_ind - 0.5


def get_binned_modProbs_near_roi(in_refPos_modProbs, roi_coord):
    if len(in_refPos_modProbs) == 0:
        return []
    vec_refPos, vec_modProb = np.vstack(in_refPos_modProbs).T
    diff_refPos = vec_refPos - roi_coord
    out_binned_modProbs = []
    for bin_ind in range(len(roi_bin_edges)-1):
        bin_start = roi_bin_edges[bin_ind]
        bin_end = roi_bin_edges[bin_ind+1]
        out_binned_modProbs.append(vec_modProb[(diff_refPos >= bin_start) * (diff_refPos < bin_end)])
    return out_binned_modProbs


def aggregate_modProb_bins(in_modProbs):
    out_modProbs = {this_mod: [[] for i in range(len(roi_bin_edges)-1)] for this_mod in mods}
    for this_mod in mods:
        for this_read_bins in in_modProbs[this_mod]:
            for i, this_bin in enumerate(this_read_bins):
                out_modProbs[this_mod][i].extend(this_bin)
    return out_modProbs


def calc_log2fc(mod_probs_1, mod_probs_2):
    out_log2fc = []
    for this_bin_1, this_bin_2 in zip(mod_probs_1, mod_probs_2):
        if (len(this_bin_1) == 0) or (len(this_bin_2) == 0):
            out_log2fc.append([])
        else:
            n1 = np.mean(np.array(this_bin_1) >= 0.5)
            n2 = np.mean(np.array(this_bin_2) >= 0.5)
            if (n1 == 0) or (n2 == 0):
                out_log2fc.append([])
            else:
                out_log2fc.append([np.log2(n1/n2)])
    return out_log2fc


thresh_IntronDepth = 5
thresh_IRratio = 0.05
thresh_intronic_percentage = 0.01

res_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/{cond}'
irf_file = os.path.join(res_dir, 'IRFinder/IRFinder-IR-nondir.txt')
bam_file = os.path.join(res_dir, 'chrALL.mAFiA.reads.bam')

gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
gtf_bedtools = pybedtools.BedTool(gtf_file)

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure4'
os.makedirs(img_out, exist_ok=True)

df_irf = pd.read_csv(irf_file, sep='\t', dtype={'Chr': str})
df_irf_thresh = df_irf[
    (df_irf['IntronDepth'] >= thresh_IntronDepth)
    * (df_irf['IRratio'] >= thresh_IRratio)
]
df_irf_clean = df_irf_thresh.loc[['known-exon' not in this_name.split('/')[-1] for this_name in df_irf_thresh['Name']], :]
df_irf_clean = df_irf_clean[df_irf_clean['Warnings'] == '-']
df_irf_clean.to_csv(os.path.join(img_out, f'df_irf_clean_{ds}_{cond}.tsv'), sep='\t', index=False)
print(f'{len(df_irf_clean)} IR junctions')

unique_tx = df_irf_clean['Name'].unique()

def get_norm_coord_mod_probs(in_refPos_modProbs, in_left_exon_range, in_right_exon_range, in_strand):
    left_exon_len = in_left_exon_range[1] - in_left_exon_range[0]
    intron_len = in_right_exon_range[0] - in_left_exon_range[1]
    right_exon_len = in_right_exon_range[1] - in_right_exon_range[0]

    mod_norm_pos_mod_probs = {}
    for this_mod in mods:
        norm_left_refPos_modProbs = [((this_refPos-in_left_exon_range[0])/left_exon_len, modProb)
                                     for this_refPos, modProb in in_refPos_modProbs[this_mod]
                                     if (this_refPos<=in_left_exon_range[1]) and (this_refPos>=in_left_exon_range[0])]
        norm_intron_refPos_modProbs = [((this_refPos-in_left_exon_range[1])/intron_len, modProb)
                                       for this_refPos, modProb in in_refPos_modProbs[this_mod]
                                       if (this_refPos<=in_right_exon_range[0]) and (this_refPos>=in_left_exon_range[1])]
        norm_right_refPos_modProbs = [((this_refPos-in_right_exon_range[0])/right_exon_len, modProb)
                                      for this_refPos, modProb in in_refPos_modProbs[this_mod]
                                      if (this_refPos<=in_right_exon_range[1]) and (this_refPos>=in_right_exon_range[0])]
        if in_strand == '+':
            mod_norm_pos_mod_probs[this_mod] = norm_left_refPos_modProbs \
                                               + [(x+1.0, p) for x, p in norm_intron_refPos_modProbs] \
                                               + [(x+2.0, p) for x, p in norm_right_refPos_modProbs]
        else:
            mod_norm_pos_mod_probs[this_mod] = [(1.0-x, p) for x, p in norm_right_refPos_modProbs] \
                                               + [(2.0-x, p) for x, p in norm_intron_refPos_modProbs] \
                                               + [(3.0-x, p) for x, p in norm_left_refPos_modProbs]
    return mod_norm_pos_mod_probs


num_bins = 5
bin_edges = np.concatenate([
    np.linspace(0, 1, num_bins+1)[:-1], np.linspace(1, 2, num_bins+1)[:-1], np.linspace(2, 3, num_bins+1)
    ])
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])


def get_binned_mod_prob_profile(in_mod_norm_pos_mod_probs):
    mod_profile = {}
    for this_mod in mods:
        this_mod_pos_probs = [this_tup for this_dict in in_mod_norm_pos_mod_probs for this_tup in this_dict[this_mod]]
        if len(this_mod_pos_probs):
            vec_pos, vec_prob = np.vstack(this_mod_pos_probs).T
            this_mod_profile = []
            for bin_i in range(len(bin_edges)-1):
                bin_start = bin_edges[bin_i]
                bin_end = bin_edges[bin_i+1]
                binned_probs = vec_prob[(vec_pos>=bin_start) * (vec_pos<bin_end)]
                if len(binned_probs):
                    this_mod_profile.append(np.mean(binned_probs >= 0.5))
                else:
                    this_mod_profile.append(0.0)
            mod_profile[this_mod] = np.nan_to_num(this_mod_profile)

    return mod_profile


intronic_reads = []
exonic_reads = []
agg_intronic_profile = []
agg_exonic_profile = []
for _, this_row in tqdm(df_irf_clean.iterrows()):
    flag_required = 0 if this_row['Strand'] == '+' else 16

    intronic_norm_pos_mod_probs = []
    exonic_norm_pos_mod_probs = []
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        chrom, irStart, irEnd = this_row[['Chr', 'Start', 'End']].values
        if this_row['Strand'] == '+':
            roiStart = irStart
            roiEnd = irEnd
        else:
            roiStart = irEnd
            roiEnd = irStart

        left_exon_bed = pybedtools.BedTool(f'{chrom} {irStart-1} {irStart}', from_string=True)
        df_left_intersect = gtf_bedtools.intersect(left_exon_bed, wa=True).to_dataframe()
        df_left_intersect_exons = df_left_intersect[
            (df_left_intersect['feature'] == 'exon')
            * (df_left_intersect['strand'] == this_row['Strand'])
            * (df_left_intersect['end'] == irStart)
        ]
        left_exon_range = [df_left_intersect_exons['start'].max(), df_left_intersect_exons['end'].values[0]]

        right_exon_bed = pybedtools.BedTool(f'{chrom} {irEnd} {irEnd+1}', from_string=True)
        df_right_intersect = gtf_bedtools.intersect(right_exon_bed, wa=True).to_dataframe()
        df_right_intersect_exons = df_right_intersect[
            (df_right_intersect['feature'] == 'exon')
            * (df_right_intersect['strand'] == this_row['Strand'])
            * (df_right_intersect['start'] == irEnd+1)
        ]
        right_exon_range = [df_right_intersect_exons['start'].values[0], df_right_intersect_exons['end'].min()]

        for this_read in bam.fetch(chrom, irStart, irEnd):
            if this_read.flag == flag_required:
                ref_to_read_pos = {ref_pos: read_pos for read_pos, ref_pos in this_read.get_aligned_pairs()}
                read_to_ref_pos = {v: k for k, v in ref_to_read_pos.items()}
                refPos_modProbs = {
                    this_mod:
                        [(read_to_ref_pos[read_pos], mod_prob / 255.0)
                         for read_pos, mod_prob in this_read.modified_bases.get(('N', 0, dict_mod_code[this_mod]), [])]
                    for this_mod in mods
                }
                intronic_read_pos = [ref_to_read_pos[this_pos] for this_pos in range(irStart, irEnd)
                                     if ref_to_read_pos.get(this_pos) is not None]
                intronic_percentage = len(intronic_read_pos) / (irEnd - irStart)
                left_exonic_read_pos = [ref_to_read_pos[this_pos] for this_pos in range(*left_exon_range)
                                        if ref_to_read_pos.get(this_pos) is not None]
                right_exonic_read_pos = [ref_to_read_pos[this_pos] for this_pos in range(*right_exon_range)
                                        if ref_to_read_pos.get(this_pos) is not None]
                if len(left_exonic_read_pos) and len(right_exonic_read_pos):
                    if intronic_percentage >= thresh_intronic_percentage:
                        intronic_reads.append(this_read)
                        intronic_norm_pos_mod_probs.append(
                            get_norm_coord_mod_probs(refPos_modProbs, left_exon_range, right_exon_range, this_row['Strand'])
                        )
                    else:
                        exonic_reads.append(this_read)
                        exonic_norm_pos_mod_probs.append(
                            get_norm_coord_mod_probs(refPos_modProbs, left_exon_range, right_exon_range, this_row['Strand'])
                        )


    agg_intronic_profile.append(get_binned_mod_prob_profile(intronic_norm_pos_mod_probs))
    agg_exonic_profile.append(get_binned_mod_prob_profile(exonic_norm_pos_mod_probs))

total_intronic_profile = {this_mod: np.mean(np.vstack([this_profile[this_mod]
                                     for this_profile in agg_intronic_profile
                                     if this_profile.get(this_mod) is not None]), axis=0)
                          for this_mod in mods}
total_exonic_profile = {this_mod: np.mean(np.vstack([this_profile[this_mod]
                                   for this_profile in agg_exonic_profile
                                   if this_profile.get(this_mod) is not None]), axis=0)
                        for this_mod in mods}

splice_type = ['spliced', 'IR']

total_profile = {
    'IR': total_intronic_profile,
    'spliced': total_exonic_profile
}

splice_colors = {
    'IR': 'r',
    'spliced': 'b'
}

ymax = np.round(
    (np.max(np.concatenate(
        [sub_val for val in total_profile.values() for sub_val in val.values()])
    ) // 0.05 + 1) * 0.05, 2
)

plt.figure(figsize=(4*cm, 4*cm))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(2, 1, mod_ind+1)
    for shift, this_splice in enumerate(splice_type):
        plt.plot(bin_centers, total_profile[this_splice][this_mod], f'{splice_colors[this_splice]}-', label=this_splice)
        plt.plot(bin_centers, total_profile[this_splice][this_mod], f'{splice_colors[this_splice]}.', markersize=1)
    plt.ylim([-0.01, ymax+0.01])
    plt.yticks(np.arange(0, ymax+0.01, 0.05))
    # plt.ylabel(rf'$\langle$$S_{{{dict_mod_display[this_mod]}}}$$\rangle$', fontsize=15)
    # plt.legend(loc='upper left')
    plt.axvspan(xmin=1, xmax=2, fc='gray', alpha=0.5)
    if mod_ind==1:
        plt.xticks([0.5, 1.5, 2.5], ['5\' exon', 'intron', '3\' exon'])
    else:
        plt.xticks([])
# plt.suptitle(f'{ds}, {cond}\nAvg. profile of {len(df_irf_clean)} IR junctions', fontsize=20)
plt.savefig(os.path.join(img_out, f'mod_profile_IR_spliced_{ds}_{cond}.{FMT}'), **fig_kwargs)

########################################################################################################################
### output bam file ####################################################################################################
########################################################################################################################
# intronic_bamfile = os.path.join(res_dir, 'IRFinder', 'intronic_reads.bam')
# exonic_bamfile = os.path.join(res_dir, 'IRFinder', 'exonic_reads.bam')
# written_intronic_read_ids = []
# written_exonic_read_ids = []
# with pysam.AlignmentFile(bam_file, 'rb') as in_bam:
#     with pysam.AlignmentFile(intronic_bamfile, 'wb', template=in_bam) as out_bam:
#         for this_read in intronic_reads:
#             if this_read.query_name not in written_intronic_read_ids:
#                 out_bam.write(this_read)
#                 written_intronic_read_ids.append(this_read.query_name)
#     with pysam.AlignmentFile(exonic_bamfile, 'wb', template=in_bam) as out_bam:
#         for this_read in exonic_reads:
#             if this_read.query_name not in written_exonic_read_ids:
#                 out_bam.write(this_read)
#                 written_exonic_read_ids.append(this_read.query_name)
#
# pysam.sort("-o", intronic_bamfile.replace('.bam', '.sorted.bam'), intronic_bamfile)
# os.rename(intronic_bamfile.replace('.bam', '.sorted.bam'), intronic_bamfile)
# pysam.index(intronic_bamfile)
# pysam.sort("-o", exonic_bamfile.replace('.bam', '.sorted.bam'), exonic_bamfile)
# os.rename(exonic_bamfile.replace('.bam', '.sorted.bam'), exonic_bamfile)
# pysam.index(exonic_bamfile)

########################################################################################################################
### compare mean delta #################################################################################################
########################################################################################################################

# mean_delta_roiStart_modProbs = aggregate_modProb_bins(mean_delta_roiStart_modProbs)
# mean_delta_roiEnd_modProbs = aggregate_modProb_bins(mean_delta_roiEnd_modProbs)
#
# mean_delta_modProbs = {
#     'start': mean_delta_roiStart_modProbs,
#     'end': mean_delta_roiEnd_modProbs
# }
#
# ylim = [-10, 10]
# plt.figure(figsize=(10, 4))
# for col_ind, roi in enumerate(['start', 'end']):
#     plt.subplot(1, 2, col_ind+1)
#     for mod_ind, this_mod in enumerate(mods):
#         for mini_col_ind, this_delta in enumerate(mean_delta_roiStart_modProbs[this_mod]):
#             this_delta = [x for x in this_delta if ~np.isinf(x)]
#             bplot = plt.boxplot(this_delta, positions=[mini_col_ind-0.2+0.4*mod_ind],
#                                 showfliers=False, patch_artist=True, whis=1.5)
#             for patch in bplot['boxes']:
#                 patch.set_facecolor(mod_colors[this_mod])
#             # for patch in bplot['means']:
#             #     patch.set_color('k')
#             for patch in bplot['medians']:
#                 patch.set_color('k')
#         if roi == 'start':
#             plt.xlim([-0.5, x_zero])
#             plt.xticks(np.arange(x_zero + 1) - 0.5, roi_bin_edges[:(x_zero_ind + 1)])
#         elif roi == 'end':
#             plt.xlim([x_zero, len(roi_bin_edges)-1.5])
#             plt.xticks(np.arange(x_zero, len(roi_bin_edges)-1), roi_bin_edges[x_zero_ind:])
#         plt.axhline(y=0, color='gray', alpha=0.5)
#         plt.ylim(ylim)
#         # plt.title(f'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
#         plt.legend(*zip(*mod_labels), loc='upper left')
#         plt.xlabel(f'Distance from IR {roi}', fontsize=12)
#         plt.ylabel(r'log2fc ($N_{IR}$ / $N_{spliced}$)', fontsize=12)
# plt.suptitle(f'{ds}, {cond}', fontsize=20)
# plt.savefig(os.path.join(img_out, f'log2gc_by_IR_distance_{ds}_{cond}.png'), bbox_inches='tight')

########################################################################################################################
### pileup #############################################################################################################
########################################################################################################################
# mod_sites = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.bed'
# os.system(
#     "python3 ~/git/mAFiA_dev/mAFiA/mAFiA_pileup.py"
#     f" --bam_file {intronic_bamfile}"
#     f" --mod_file {mod_sites}"
#     " --min_coverage 1"
#     f" --out_dir {os.path.join(res_dir, 'IRFinder')}"
#     f" --out_filename {intronic_bamfile.replace('.bam', '.bed')}"
# )

########################################################################################################################
### avg mod per read ###################################################################################################
########################################################################################################################
# intronic_mod_level = calc_mod_level_per_read(intronic_reads)
# exonic_mod_level = calc_mod_level_per_read(exonic_reads)
#
# num_bins = 50
# bin_edges = np.linspace(0, 1, num_bins+1)
# bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
# ylim = [0, 0.08]
#
# plt.figure(figsize=(10, 4))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     intronic_vals = [this_val[this_mod] for this_val in intronic_mod_level.values()]
#     hist_intron, _ = np.histogram(intronic_vals, bins=bin_edges)
#     norm_hist_intron = hist_intron / np.sum(hist_intron)
#     exonic_vals = [this_val[this_mod] for this_val in exonic_mod_level.values()]
#     hist_exon, _ = np.histogram(exonic_vals, bins=bin_edges)
#     norm_hist_exon = hist_exon / np.sum(hist_exon)
#     plt.plot(bin_centers, norm_hist_intron, c='r', label=f'intronic ({len(intronic_vals)})')
#     plt.plot(bin_centers, norm_hist_exon, c='b', label=f'non-intronic ({len(exonic_vals)})')
#     plt.legend()
#     plt.ylim(ylim)
#     # plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $', fontsize=12)
#     # plt.ylabel('Density', fontsize=12)
#     # plt.title(f'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
# # plt.suptitle(f'{ds}, {cond}', fontsize=20)
# plt.savefig(os.path.join(img_out, f'hist_IR_versus_spliced_{ds}_{cond}.{FMT}'), **fig_kwargs)