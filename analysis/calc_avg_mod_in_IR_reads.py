import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from tqdm import tqdm
import pysam
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pybedtools


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

# roi_bin_edges = [-100, -50, 0, 50, 100]
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
            n1 = np.sum(np.array(this_bin_1) >= 0.5)
            n2 = np.sum(np.array(this_bin_2) >= 0.5)
            if (n1 == 0) or (n2 == 0):
                out_log2fc.append([])
            else:
                out_log2fc.append([np.log2(n1/n2)])
    return out_log2fc


thresh_IntronDepth = 10
thresh_IRratio = 0.1
thresh_intronic_percentage = 0.1

ds = 'HFpEF'
cond = 'ctrl_merged'
res_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/{cond}'
irf_file = os.path.join(res_dir, 'IRFinder/IRFinder-IR-nondir.txt')
bam_file = os.path.join(res_dir, 'chrALL.mAFiA.reads.bam')

gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
gtf_bedtools = pybedtools.BedTool(gtf_file)

img_out = '/home/adrian/img_out/IR_reads_avg_mod'
os.makedirs(img_out, exist_ok=True)

df_irf = pd.read_csv(irf_file, sep='\t', dtype={'Chr': str})
df_irf_thresh = df_irf[
    (df_irf['IntronDepth'] >= thresh_IntronDepth)
    * (df_irf['IRratio'] >= thresh_IRratio)
    * (df_irf['IRratio'] < 1.0)
]
# df_irf_clean = df_irf_thresh
df_irf_clean = df_irf_thresh.loc[[this_name.split('/')[-1] == 'clean' for this_name in df_irf_thresh['Name']], :]
# df_irf_clean = df_irf_clean[df_irf_clean['Warnings'] == '-']

intronic_reads = []
exonic_reads = []
mean_delta_roiStart_modProbs = {this_mod: [] for this_mod in mods}
mean_delta_roiEnd_modProbs = {this_mod: [] for this_mod in mods}
for _, this_row in tqdm(df_irf_clean.iterrows()):
# for this_name in tqdm(df_irf_clean['Name'].unique()):
#     this_row = df_irf_clean[df_irf_clean['Name'] == this_name].iloc[0]
    intronic_roiStart_modProbs = {this_mod: [] for this_mod in mods}
    intronic_roiEnd_modProbs = {this_mod: [] for this_mod in mods}
    exonic_roiStart_modProbs = {this_mod: [] for this_mod in mods}
    exonic_roiEnd_modProbs = {this_mod: [] for this_mod in mods}
    flag_required = 0 if this_row['Strand'] == '+' else 16
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        chrom, irStart, irEnd = this_row[['Chr', 'Start', 'End']].values
        if this_row['Strand'] == '+':
            roiStart = irStart
            roiEnd = irEnd
        else:
            roiStart = irEnd
            roiEnd = irStart

        # this_row_bed = pybedtools.BedTool(f'{chrom} {irStart} {irEnd}', from_string=True)
        # intersect = gtf_bedtools.intersect(this_row_bed, wa=True)

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
                intronic_read_pos = [ref_to_read_pos.get(this_pos) for this_pos in range(irStart, irEnd)]
                intronic_percentage = len([x for x in intronic_read_pos if x is not None]) / len(intronic_read_pos)
                if intronic_percentage >= thresh_intronic_percentage:
                    intronic_reads.append(this_read)
                    for this_mod in mods:
                        intronic_roiStart_modProbs[this_mod].append(
                            get_binned_modProbs_near_roi(refPos_modProbs[this_mod], roiStart))
                        intronic_roiEnd_modProbs[this_mod].append(
                            get_binned_modProbs_near_roi(refPos_modProbs[this_mod], roiEnd))
                else:
                    exonic_reads.append(this_read)
                    for this_mod in mods:
                        exonic_roiStart_modProbs[this_mod].append(
                            get_binned_modProbs_near_roi(refPos_modProbs[this_mod], roiStart))
                        exonic_roiEnd_modProbs[this_mod].append(
                            get_binned_modProbs_near_roi(refPos_modProbs[this_mod], roiEnd))

    intronic_roiStart_modProbs = aggregate_modProb_bins(intronic_roiStart_modProbs)
    exonic_roiStart_modProbs = aggregate_modProb_bins(exonic_roiStart_modProbs)
    intronic_roiEnd_modProbs = aggregate_modProb_bins(intronic_roiEnd_modProbs)
    exonic_roiEnd_modProbs = aggregate_modProb_bins(exonic_roiEnd_modProbs)

    for this_mod in mods:
        mean_delta_roiStart_modProbs[this_mod].append(
            calc_log2fc(intronic_roiStart_modProbs[this_mod], exonic_roiStart_modProbs[this_mod])
        )
        mean_delta_roiEnd_modProbs[this_mod].append(
            calc_log2fc(intronic_roiEnd_modProbs[this_mod], exonic_roiEnd_modProbs[this_mod])
        )

mean_delta_roiStart_modProbs = aggregate_modProb_bins(mean_delta_roiStart_modProbs)
mean_delta_roiEnd_modProbs = aggregate_modProb_bins(mean_delta_roiEnd_modProbs)

mean_delta_modProbs = {
    'start': mean_delta_roiStart_modProbs,
    'end': mean_delta_roiEnd_modProbs
}

ylim = [-10, 10]
plt.figure(figsize=(10, 4))
for col_ind, roi in enumerate(['start', 'end']):
    plt.subplot(1, 2, col_ind+1)
    for mod_ind, this_mod in enumerate(mods):
        for mini_col_ind, this_delta in enumerate(mean_delta_roiStart_modProbs[this_mod]):
            this_delta = [x for x in this_delta if ~np.isinf(x)]
            bplot = plt.boxplot(this_delta, positions=[mini_col_ind-0.2+0.4*mod_ind],
                                showfliers=False, patch_artist=True, whis=1.5)
            for patch in bplot['boxes']:
                patch.set_facecolor(mod_colors[this_mod])
            # for patch in bplot['means']:
            #     patch.set_color('k')
            for patch in bplot['medians']:
                patch.set_color('k')
        if roi == 'start':
            plt.xlim([-0.5, x_zero])
            plt.xticks(np.arange(x_zero + 1) - 0.5, roi_bin_edges[:(x_zero_ind + 1)])
        elif roi == 'end':
            plt.xlim([x_zero, len(roi_bin_edges)-1.5])
            plt.xticks(np.arange(x_zero, len(roi_bin_edges)-1), roi_bin_edges[x_zero_ind:])
        plt.axhline(y=0, color='gray', alpha=0.5)
        plt.ylim(ylim)
        # plt.title(f'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
        plt.legend(*zip(*mod_labels), loc='upper left')
        plt.xlabel(f'Distance from IR {roi}', fontsize=12)
        plt.ylabel(r'log2fc ($N_{IR}$ / $N_{non-IR}$)', fontsize=12)
plt.suptitle(f'{ds}, {cond}', fontsize=20)
plt.savefig(os.path.join(img_out, f'log2gc_by_IR_distance_{ds}_{cond}.png'), bbox_inches='tight')

########################################################################################################################
### output bam file ####################################################################################################
########################################################################################################################
intronic_bamfile = os.path.join(res_dir, 'IRFinder', 'intronic_reads.bam')
exonic_bamfile = os.path.join(res_dir, 'IRFinder', 'exonic_reads.bam')
with pysam.AlignmentFile(bam_file, 'rb') as in_bam:
    with pysam.AlignmentFile(intronic_bamfile, 'wb', template=in_bam) as out_bam:
        for this_read in intronic_reads:
            out_bam.write(this_read)
    with pysam.AlignmentFile(exonic_bamfile, 'wb', template=in_bam) as out_bam:
        for this_read in exonic_reads:
            out_bam.write(this_read)

pysam.sort("-o", intronic_bamfile.replace('.bam', '.sorted.bam'), intronic_bamfile)
os.rename(intronic_bamfile.replace('.bam', '.sorted.bam'), intronic_bamfile)
pysam.index(intronic_bamfile)
pysam.sort("-o", exonic_bamfile.replace('.bam', '.sorted.bam'), exonic_bamfile)
os.rename(exonic_bamfile.replace('.bam', '.sorted.bam'), exonic_bamfile)
pysam.index(exonic_bamfile)

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
intronic_mod_level = calc_mod_level_per_read(intronic_reads)
exonic_mod_level = calc_mod_level_per_read(exonic_reads)

num_bins = 50
bin_edges = np.linspace(0, 1, num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
ylim = [0, 0.08]

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    intronic_vals = [this_val[this_mod] for this_val in intronic_mod_level.values()]
    hist_intron, _ = np.histogram(intronic_vals, bins=bin_edges)
    norm_hist_intron = hist_intron / np.sum(hist_intron)
    exonic_vals = [this_val[this_mod] for this_val in exonic_mod_level.values()]
    hist_exon, _ = np.histogram(exonic_vals, bins=bin_edges)
    norm_hist_exon = hist_exon / np.sum(hist_exon)
    plt.plot(bin_centers, norm_hist_intron, c='r', label=f'intronic ({len(intronic_vals)})')
    plt.plot(bin_centers, norm_hist_exon, c='b', label=f'non-intronic ({len(exonic_vals)})')
    plt.legend()
    plt.ylim(ylim)
    plt.xlabel(rf'$\langle P({{{dict_mod_display[this_mod]}}}) \rangle $', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title(f'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
plt.suptitle(f'{ds}, {cond}', fontsize=20)
# plt.savefig(os.path.join(img_out, f'hist_{ds}_{cond}.png'), bbox_inches='tight')