import os
import pandas as pd
from scipy.spatial.distance import pdist
import pysam
import numpy as np
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 6
mpl.rcParams['legend.fontsize'] = 3
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
######################################################################
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm

# img_out = '/home/adrian/img_out/single_read_cross_talk'
img_out = '/home/achan/img_out/single_read_cross_talk'
os.makedirs(img_out, exist_ok=True)
# img_out = '/home/adrian/img_out/manuscript_bioinformatics_application_note'
# os.makedirs(img_out, exist_ok=True)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}


dist_bin_max = 1000
dist_num_bins = 200
dist_bin_edges = np.linspace(-dist_bin_max, dist_bin_max, dist_num_bins+1)
dist_bin_centers = 0.5 * (dist_bin_edges[1:] + dist_bin_edges[:-1])

THRESH_PROB = 0.5

def get_hist_dist_from_loc_prob(loc_prob_1, loc_prob_2=None, thresh_prob=0.5):
    thresh_cross_dist = []
    vec_loc_1 = np.array([loc_prob[0] for loc_prob in loc_prob_1])

    if loc_prob_2 is None:
        if len(vec_loc_1):
            dist_all = pdist(vec_loc_1[:, np.newaxis])
            vec_loc_thresh_1 = np.array([loc_prob[0] for loc_prob in loc_prob_1 if loc_prob[1] >= thresh_prob])
            if len(vec_loc_thresh_1):
                dist_thresh = pdist(vec_loc_thresh_1[:, np.newaxis])
            else:
                dist_thresh = []
        else:
            dist_all = []
            dist_thresh = []

    else:
        vec_loc_2 = np.array([loc_prob[0] for loc_prob in loc_prob_2])

        if len(vec_loc_1) and len(vec_loc_2):
            dist_all = (vec_loc_1[:, np.newaxis] - vec_loc_2[np.newaxis, :]).flatten()
            vec_loc_thresh_1 = np.array([loc_prob[0] for loc_prob in loc_prob_1 if loc_prob[1] >= thresh_prob])
            vec_loc_thresh_2 = np.array([loc_prob[0] for loc_prob in loc_prob_2 if loc_prob[1] >= thresh_prob])
            if len(vec_loc_thresh_1) and len(vec_loc_thresh_2):
                dist_thresh = (vec_loc_thresh_1[:, np.newaxis] - vec_loc_thresh_2[np.newaxis, :]).flatten()
            else:
                dist_thresh = []
        else:
            dist_all = []
            dist_thresh = []

    dist_all = [this_dist for this_dist in dist_all
                if (this_dist >= -dist_bin_max) and (this_dist <= dist_bin_max)]
    dist_thresh = [this_dist for this_dist in dist_thresh
                   if (this_dist >= -dist_bin_max) and (this_dist <= dist_bin_max)]

    if len(dist_all):
        hist_dist_all, _ = np.histogram(dist_all, bins=dist_bin_edges)
    else:
        hist_dist_all = None
    if len(dist_thresh):
        hist_dist_thresh, _ = np.histogram(dist_thresh, bins=dist_bin_edges)
    else:
        hist_dist_thresh = None

    return hist_dist_all, hist_dist_thresh


def get_read_mod_distances(in_read, thresh_prob):
    read_to_ref_loc = {read_loc: ref_loc for read_loc, ref_loc in in_read.get_aligned_pairs(matches_only=True)}
    mod_loc_prob = {}
    for this_mod, this_tag in mod_tags.items():
        raw_loc_prob = in_read.modified_bases.get(this_tag, [])
        if len(raw_loc_prob):
            mod_loc_prob[this_mod] = \
                [(read_to_ref_loc[loc], prob / 255.0) for loc, prob in raw_loc_prob if loc in read_to_ref_loc.keys()]

    mod_dist_hist = {}
    for this_mod in mod_tags.keys():
        if this_mod in mod_loc_prob.keys():
            mod_dist_hist[this_mod] = get_hist_dist_from_loc_prob(mod_loc_prob[this_mod], thresh_prob)
        else:
            mod_dist_hist[this_mod] = (None, None)

    if ('m6A' in mod_loc_prob.keys()) and ('psi' in mod_loc_prob.keys()):
        mod_dist_hist['cross'] = get_hist_dist_from_loc_prob(mod_loc_prob['m6A'], mod_loc_prob['psi'])
    else:
        mod_dist_hist['cross'] = (None, None)

    return mod_dist_hist

########################################################################################################################
### R002 ###############################################################################################################
########################################################################################################################
mod_tags = {
    'm6A': ('N', 0, 21891),
    'psi': ('N', 0, 17802)
}

# base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'
base_dir = '/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'


ds = 'WT'
bam_file = os.path.join(base_dir, 'HEK293/WT_P2/chrALL.mAFiA.reads.bam')

# ds = 'M3KO'
# bam_file = os.path.join(base_dir, 'HEK293T_Mettl3_KO/merged/chrALL.mAFiA.reads.bam')

# ds = 'M3KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siMETTL3_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siTRUB1_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1OE'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_OE/merged/chrALL.mAFiA.reads.bam')

########################################################################################################################
### R004 ###############################################################################################################
########################################################################################################################
# mod_tags = {
#     'm6A': ('A', 0, 'a'),
#     'psi': ('T', 0, 17802)
# }
#
# base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'
#
# ds = 'HEK293_ctrl_R004'
# bam_file = os.path.join(base_dir, 'HEK293_ctrl_RTA/calls_2025-02-26_T06-44-51.bam')

# ds = 'HEK293_TRUB1_kd'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_kd_RTA/calls_2025-02-26_T06-43-59.bam')

########################################################################################################################

outfile_name = os.path.join(img_out, f'read_level_mod_distance_{ds}.{FMT}')
# outfile_name = os.path.join(img_out, f'figureS4.{FMT}')

# gene_bed = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/gene.ensembl_havana.GRCh38.102.bed'
# df_gene = pd.read_csv(gene_bed, sep='\t')

single_read_mod_distances = []
with pysam.AlignmentFile(bam_file, 'rb', check_sq=False) as bam:
    for this_read in tqdm(bam.fetch(until_eof=True)):
        if this_read.modified_bases is not None:
            single_read_mod_distances.append(get_read_mod_distances(this_read), THRESH_PROB)

dict_mod_display['cross'] = '\psi$$\mapsto$$m^6A'

for mod_type in ['m6A', 'psi', 'cross']:
    hist_dist_all = [this_read_mod_distances[mod_type][0] for this_read_mod_distances in single_read_mod_distances
                     if this_read_mod_distances[mod_type][0] is not None]
    avg_hist_dist_all = np.mean(np.vstack(hist_dist_all), axis=0)
    norm_hist_dist_all = avg_hist_dist_all / np.sum(avg_hist_dist_all)
    cdf_from_center_all = [np.sum(norm_hist_dist_all[(int(dist_num_bins/2)-i-1):(int(dist_num_bins/2)+i+1)]) for i in range(int(dist_num_bins/2))]
    quartile_all = dist_bin_centers[np.where(np.array(cdf_from_center_all) >= 0.50)[0][0] + int(dist_num_bins / 2)]

    hist_dist_thresh = [this_read_mod_distances[mod_type][1] for this_read_mod_distances in single_read_mod_distances
                        if this_read_mod_distances[mod_type][1] is not None]
    avg_hist_dist_thresh = np.mean(np.vstack(hist_dist_thresh), axis=0)
    norm_hist_dist_thresh = avg_hist_dist_thresh / np.sum(avg_hist_dist_thresh)
    cdf_from_center_thresh = [np.sum(norm_hist_dist_thresh[(int(dist_num_bins/2)-i-1):(int(dist_num_bins/2)+i+1)]) for i in range(int(dist_num_bins/2))]
    quartile_thresh = dist_bin_centers[np.where(np.array(cdf_from_center_thresh) >= 0.50)[0][0] + int(dist_num_bins / 2)]

    plt.figure(figsize=(4 * cm, 4 * cm))
    plt.plot(dist_bin_centers, norm_hist_dist_all, c='b', label=f"${dict_mod_display[mod_type]}$ all\n$d_{{{50}}}$={int(quartile_all)}nts")
    plt.plot(dist_bin_centers, norm_hist_dist_thresh, c='r',
             label=f"${dict_mod_display[mod_type]}$ P$\geq${THRESH_PROB}\n$d_{{{50}}}$={int(quartile_thresh)}nts")
    if mod_type == 'cross':
        plt.xticks(np.linspace(-dist_bin_max, dist_bin_max, 5))
        plt.xlim([-dist_bin_max, dist_bin_max])
    else:
        plt.xticks(np.linspace(0, dist_bin_max, 5))
        plt.xlim([0, dist_bin_max])
    plt.xlabel('Distance (nt)')
    plt.ylabel('Probability')
    plt.legend(loc='upper right')
    plt.savefig(os.path.join(img_out, f'norm_hist_dist_{mod_type}.{FMT}'), **fig_kwargs)

# hist_dist_psi_all = [this_read_mod_distances['psi'][0] for this_read_mod_distances in single_read_mod_distances
#                      if this_read_mod_distances['psi'][0] is not None]
# hist_dist_psi_thresh = [this_read_mod_distances['psi'][1] for this_read_mod_distances in single_read_mod_distances
#                         if this_read_mod_distances['psi'][1] is not None]
# hist_dist_cross_all = [this_read_mod_distances['cross'][0] for this_read_mod_distances in single_read_mod_distances
#                      if this_read_mod_distances['cross'][0] is not None]
# hist_dist_cross_thresh = [this_read_mod_distances['cross'][1] for this_read_mod_distances in single_read_mod_distances
#                         if this_read_mod_distances['cross'][1] is not None]
#
# plt.figure(figsize=(4*cm, 4*cm))
# plt.plot(dist_bin_centers, norm_hist_dist_psi_all, c='b', label=f"${dict_mod_display['psi']}$ all")
# plt.plot(dist_bin_centers, norm_hist_dist_psi_thresh, c='r', label=f"${dict_mod_display['psi']}$ P$\geq${THRESH_PROB}")
# plt.xticks(np.linspace(0, dist_bin_max, 5))
# plt.xlim([0, dist_bin_max])
# plt.xlabel('Distance (nt)')
# plt.ylabel('Probability')
# plt.legend(loc='upper right')
# plt.savefig(os.path.join(img_out, f'norm_hist_dist_m6A.{FMT}'), **fig_kwargs)











flierprops = dict(marker='o', markerfacecolor='none', markersize=2, markeredgecolor='gray',
                  alpha=0.5, rasterized=True)

xy_ticks = np.int32(bin_edges * 100)

plt.figure(figsize=(8*cm, 7*cm))
plt.subplot(2, 2, 1)
plt.boxplot(binned_psi, flierprops=flierprops)
# plt.violinplot(binned_psi, quantiles=[[0.5]]*len(binned_psi))
# plt.plot(np.arange(1, len(bin_edges)), top_psi, c='r', ls='-')
# plt.plot(np.arange(1, len(bin_edges)), top_psi, 'r+', markersize=2, label=f'N$\geq${thresh_top_reads}')
# bin_sizes = [len(this_bin) for this_bin in binned_psi]
# for bin_ind, this_bin_size in enumerate(bin_sizes):
#     plt.text(bin_ind+0.5, 1.05, this_bin_size)
# plt.legend(loc='upper right')
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"N(${dict_mod_display['m6A']}$)")
plt.ylabel(f"N(${dict_mod_display['psi']}$)")
plt.subplot(2, 2, 2)
plt.boxplot(binned_m6A, flierprops=flierprops)
# plt.violinplot(binned_m6A, quantiles=[[0.5]]*len(binned_m6A))
# plt.plot(np.arange(1, len(bin_edges)), top_m6A, c='r', ls='-')
# plt.plot(np.arange(1, len(bin_edges)), top_m6A, 'r+', markersize=2, label=f'N$\geq${thresh_top_reads}')
# plt.legend(loc='upper right')
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"N(${dict_mod_display['psi']}$)")
plt.ylabel(f"N(${dict_mod_display['m6A']}$)")
# plt.savefig(os.path.join(img_out, f'boxplot_mean_occupancy_per_read_top{num_top_locs}_{ds}.{FMT}'), **fig_kwargs)
# plt.suptitle(f'{ds}\nMin. {thresh_min_locs} locs per read')

# thresh_top_reads = 0.75
# top_psi = [np.mean(this_bin[this_bin >= thresh_top_reads]) for this_bin in binned_psi]
# top_m6A = [np.mean(this_bin[this_bin >= thresh_top_reads]) for this_bin in binned_m6A]
top_reads = 100
label = f'Top {top_reads}'
top_psi = [np.mean(np.sort(this_bin)[-top_reads:]) for this_bin in binned_psi]
top_m6A = [np.mean(np.sort(this_bin)[-top_reads:]) for this_bin in binned_m6A]

# ylim = [0.79, 0.85]
# yticks = np.linspace(*ylim, 3)
top_color = 'k'
# plt.figure(figsize=(8*cm, 4*cm))
plt.subplot(2, 2, 3)
plt.plot(np.arange(1, len(bin_edges)), top_psi, c=top_color, ls='-', label=label)
plt.plot(np.arange(1, len(bin_edges)), top_psi, f'{top_color}o', markersize=2)
plt.legend(loc='lower left')
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.ylim([-0.01, 1.05])
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"N(${dict_mod_display['m6A']}$)")
plt.ylabel(f"N(${dict_mod_display['psi']}$)")
# plt.ylim(ylim)
plt.subplot(2, 2, 4)
plt.plot(np.arange(1, len(bin_edges)), top_m6A, c=top_color, ls='-', label=label)
plt.plot(np.arange(1, len(bin_edges)), top_m6A, f'{top_color}o', markersize=2)
plt.legend(loc='lower left')
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.ylim([-0.01, 1.05])
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"N(${dict_mod_display['psi']}$)")
plt.ylabel(f"N(${dict_mod_display['m6A']}$)")
# plt.ylim(ylim)
# plt.suptitle(f'N$\geq${int(thresh_top_reads*100)}%')
plt.tight_layout()
# plt.savefig(os.path.join(img_out, f'boxplot_mean_occupancy_per_read_{ds}_above{thresh_top_reads}.{FMT}'), **fig_kwargs)
plt.savefig(outfile_name, **fig_kwargs)
