import os
import pandas as pd
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
mpl.rcParams['legend.fontsize'] = 4
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
#######################################################################
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm


dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresh_valid_reads = 1000
num_top_locs = 5

def get_mean_logit(in_probs):
    if len(in_probs) == 0:
        return np.nan
    top_probs = np.sort(in_probs)[-num_top_locs:]
    rescaled_probs = np.clip(np.array(top_probs) / 255.0, a_max=0.999, a_min=0.001)
    logits = np.log2(rescaled_probs / (1-rescaled_probs))
    return np.mean(logits)


def get_mean_logit_mod_level(in_read):
    mod_mean_logit = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_probs = [this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]
        mod_mean_logit[this_mod] = get_mean_logit(this_mod_probs)
    return mod_mean_logit

thresh_min_locs = 10
def get_mod_mean_occupancy(in_read, min_locs=thresh_min_locs):
    mod_mean_occupancy = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_probs = np.array([this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]) / 255.0
        if len(this_mod_probs) >= min_locs:
            mod_mean_occupancy[this_mod] = np.mean(this_mod_probs >= 0.5)
        else:
            mod_mean_occupancy[this_mod] = np.nan
    return mod_mean_occupancy

########################################################################################################################
### R002 ###############################################################################################################
########################################################################################################################
mod_tags = {
    'm6A': ('N', 0, 21891),
    'psi': ('N', 0, 17802)
}

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'

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

# base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'

# ds = 'HEK293_ctrl_R004'
# bam_file = os.path.join(base_dir, 'HEK293_ctrl_RTA/calls_2025-02-26_T06-44-51.bam')

# ds = 'HEK293_TRUB1_kd'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_kd_RTA/calls_2025-02-26_T06-43-59.bam')

########################################################################################################################
img_out = '/home/adrian/img_out/single_read_cross_talk'
os.makedirs(img_out, exist_ok=True)

# gene_bed = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/gene.ensembl_havana.GRCh38.102.bed'
# df_gene = pd.read_csv(gene_bed, sep='\t')

single_read_mean_occupancy = []
with pysam.AlignmentFile(bam_file, 'rb', check_sq=False) as bam:
    for this_read in tqdm(bam.fetch(until_eof=True)):
        if this_read.modified_bases is not None:
            # single_read_mean_logit.append(get_mean_logit_mod_level(this_read))
            single_read_mean_occupancy.append(get_mod_mean_occupancy(this_read))

vec_m6A, vec_psi = np.vstack([
    (this_read_mean_occupancy['m6A'],  this_read_mean_occupancy['psi'])
    for this_read_mean_occupancy in single_read_mean_occupancy
    if ~np.isnan(this_read_mean_occupancy['m6A']) and ~np.isnan(this_read_mean_occupancy['psi'])
]).T
num_valid_reads = len(vec_m6A)

binned_psi = []
binned_m6A = []
bin_edges = np.round(np.linspace(0, 1, 6), 1)
for bin_i in range(len(bin_edges)-1):
    bin_start = bin_edges[bin_i]
    bin_end = bin_edges[bin_i+1]
    if bin_end == 1.0:
        bin_end += 0.001
    mask_m6A = (vec_m6A >= bin_start) * (vec_m6A < bin_end)
    binned_psi.append(vec_psi[mask_m6A])
    mask_psi = (vec_psi >= bin_start) * (vec_psi < bin_end)
    binned_m6A.append(vec_m6A[mask_psi])

top_reads = 100

flierprops = dict(marker='o', markerfacecolor='none', markersize=2, markeredgecolor='black')

xy_ticks = np.int32(bin_edges * 100)

plt.figure(figsize=(8*cm, 5*cm))
plt.subplot(1, 2, 1)
plt.boxplot(binned_psi, flierprops=flierprops)
# plt.violinplot(binned_psi, quantiles=[[0.5]]*len(binned_psi))
top_psi = [np.mean(np.sort(this_bin)[-top_reads:]) for this_bin in binned_psi]
plt.plot(np.arange(1, len(bin_edges)), top_psi, c='r', ls='-')
plt.plot(np.arange(1, len(bin_edges)), top_psi, 'r+', markersize=2, label=f'Top{top_reads} mean')
# bin_sizes = [len(this_bin) for this_bin in binned_psi]
# for bin_ind, this_bin_size in enumerate(bin_sizes):
#     plt.text(bin_ind+0.5, 1.05, this_bin_size)
plt.legend(loc='upper right')
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"N(${dict_mod_display['m6A']}$)")
plt.ylabel(f"N(${dict_mod_display['psi']}$)")
plt.subplot(1, 2, 2)
plt.boxplot(binned_m6A, flierprops=flierprops)
# plt.violinplot(binned_m6A, quantiles=[[0.5]]*len(binned_m6A))
top_m6A = [np.mean(np.sort(this_bin)[-100:]) for this_bin in binned_m6A]
plt.plot(np.arange(1, len(bin_edges)), top_m6A, c='r', ls='-')
plt.plot(np.arange(1, len(bin_edges)), top_m6A, 'r+', markersize=2, label=f'Top{top_reads} mean')
plt.legend(loc='upper right')
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"N(${dict_mod_display['psi']}$)")
plt.ylabel(f"N(${dict_mod_display['m6A']}$)")
# plt.savefig(os.path.join(img_out, f'boxplot_mean_occupancy_per_read_top{num_top_locs}_{ds}.{FMT}'), **fig_kwargs)
plt.suptitle(f'{ds}\nMin. {thresh_min_locs} locs per read')
plt.tight_layout()
plt.savefig(os.path.join(img_out, f'boxplot_mean_occupancy_per_read_{ds}.{FMT}'), **fig_kwargs)


# plt.figure(figsize=(5*cm, 5*cm))
# plt.scatter(vec_m6A, vec_psi)
# plt.savefig(os.path.join(img_out, f'scatter_mean_occupancy_per_read_top{num_top_locs}_{ds}.{FMT}'), **fig_kwargs)


# num_hh = np.sum((vec_m6A > 0) * (vec_psi > 0))
# num_hl = np.sum((vec_m6A > 0) * (vec_psi < 0))
# num_lh = np.sum((vec_m6A < 0) * (vec_psi > 0))
# num_ll = np.sum((vec_m6A < 0) * (vec_psi < 0))
# num_total = num_hh + num_hl + num_lh + num_ll
#
# perc_hh = round(num_hh / num_total * 100, 2)
# perc_hl = round(num_hl / num_total * 100, 2)
# perc_lh = round(num_lh / num_total * 100, 2)
# perc_ll = round(num_ll / num_total * 100, 2)
#
# xy_max = 10
# num_bins = 80
# mat_z, edges_x, edges_y = np.histogram2d(vec_m6A, vec_psi, bins=num_bins, range=[[-xy_max, xy_max], [-xy_max, xy_max]])
# np.savez(os.path.join(img_out, f'mean_logit_S_per_read_top{num_top_locs}_{ds}.npz'),
#          mat_z=mat_z, edges_x=edges_x, edges_y=edges_y)
#
#
# plt.figure(figsize=(4*cm, 4*cm))
# plt.axvline(x=0, c='red', ls='--')
# plt.axhline(y=0, c='red', ls='--')
# plt.imshow(np.log10(mat_z+1), extent=[-xy_max, xy_max, -xy_max, xy_max], origin='lower', vmax=2.5)
# ax = plt.gca()
# plt.text(0.1, 0.1, f'{perc_ll}%', c='r', ha='left', va='bottom', transform=ax.transAxes)
# plt.text(0.1, 0.9, f'{perc_lh}%', c='r', ha='left', va='top', transform=ax.transAxes)
# plt.text(0.9, 0.1, f'{perc_hl}%', c='r', ha='right', va='bottom', transform=ax.transAxes)
# plt.text(0.9, 0.9, f'{perc_hh}%', c='r', ha='right', va='top', transform=ax.transAxes)
# plt.xlim([-xy_max, xy_max])
# plt.ylim([-xy_max, xy_max])
# plt.xticks(np.linspace(-xy_max, xy_max, 5))
# plt.yticks(np.linspace(-xy_max, xy_max, 5))
# plt.xlabel(rf"$\langle$logit $p({{{dict_mod_display['m6A']}}})$$\rangle$ per read")
# plt.ylabel(rf"$\langle$logit $p({{{dict_mod_display['psi']}}})$$\rangle$ per read")
# plt.title(f'{ds}\n{num_valid_reads} valid reads\nTop {num_top_locs} locs')
# plt.savefig(os.path.join(img_out, f'mean_logit_S_per_read_top{num_top_locs}_{ds}.{FMT}'), **fig_kwargs)
# plt.close('all')

# with pysam.AlignmentFile(bam_file, 'rb') as bam:
#     for _, this_row in tqdm(df_gene.iterrows()):
#         this_chrom, this_chromStart, this_chromEnd, this_strand = this_row[
#             ['chrom', 'chromStart', 'chromEnd', 'strand']
#         ]
#         this_gene_id, this_gene_name = this_row[['gene_id', 'gene_name']]
#
#         flag_required = 0 if this_strand == '+' else 16
#
#         total_counts = 0
#         valid_reads = []
#         for this_read in bam.fetch(contig=this_chrom, start=this_chromStart, stop=this_chromEnd):
#             total_counts += 1
#             if this_read.flag == flag_required:
#                 valid_reads.append(this_read)
#
#         if len(valid_reads) >= thresh_valid_reads:
#             mean_logit = []
#             for this_valid_read in valid_reads:
#                 mean_logit.append(get_mean_logit_mod_level(this_valid_read))


            # vec_m6A, vec_psi = np.vstack([
            #     (this_read_mean_logit['m6A'],  this_read_mean_logit['psi'])
            #     for this_read_mean_logit in mean_logit
            # ]).T

            # xy_max = 8
            # plt.figure(figsize=(4, 4))
            # plt.axvline(x=0, c='gray', ls='--')
            # plt.axhline(y=0, c='gray', ls='--')
            # plt.scatter(vec_m6A, vec_psi, s=1)
            # plt.xlim([-xy_max, xy_max])
            # plt.ylim([-xy_max, xy_max])
            # plt.xticks(np.linspace(-xy_max, xy_max, 5))
            # plt.yticks(np.linspace(-xy_max, xy_max, 5))
            # plt.xlabel(rf"$\langle$logit ${{{dict_mod_display['m6A']}}}$$\rangle$ per read")
            # plt.ylabel(rf"$\langle$logit ${{{dict_mod_display['psi']}}}$$\rangle$ per read")
            # plt.title(f'{this_gene_id}\n{this_gene_name}')
            # plt.savefig(os.path.join(img_out, f'{this_gene_id}_{this_gene_name}.png'), bbox_inches='tight')
            # plt.close('all')

### compare 2 ds ###
# ds0 = 'WT'
# ds1 = 'TRUB1OE'
#
# data0 = np.load(os.path.join(img_out, f'mean_logit_S_per_read_top{num_top_locs}_{ds0}.npz'))
# data1 = np.load(os.path.join(img_out, f'mean_logit_S_per_read_top{num_top_locs}_{ds1}.npz'))
#
# plt.figure(figsize=(8*cm, 4*cm))
#
# plt.subplot(1, 2, 1)
# bin_edges = data0['edges_x']
# bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
# pdf0 = data0['mat_z'].sum(axis=1)
# pdf0 = pdf0 / pdf0.sum()
# pdf1 = data1['mat_z'].sum(axis=1)
# pdf1 = pdf1 / pdf1.sum()
# plt.plot(bin_centers, pdf0, label=ds0)
# plt.plot(bin_centers, pdf1, label=ds1)
# plt.title('m6A')
# plt.legend()
#
# plt.subplot(1, 2, 2)
# bin_edges = data0['edges_y']
# bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
# pdf0 = data0['mat_z'].sum(axis=0)
# pdf0 = pdf0 / pdf0.sum()
# pdf1 = data1['mat_z'].sum(axis=0)
# pdf1 = pdf1 / pdf1.sum()
# plt.plot(bin_centers, pdf0, label=ds0)
# plt.plot(bin_centers, pdf1, label=ds1)
# plt.title('psi')
# plt.legend()

