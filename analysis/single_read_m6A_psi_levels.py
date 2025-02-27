import os
import pandas as pd
import pysam
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm


mod_tags = {
    'm6A': ('N', 0, 21891),
    'psi': ('N', 0, 17802)
}
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresh_valid_reads = 1000

def get_mean_logit(in_probs):
    if len(in_probs) == 0:
        return np.nan
    rescaled_probs = np.clip(np.array(in_probs) / 255.0, a_max=0.999, a_min=0.001)
    logits = np.log2(rescaled_probs / (1-rescaled_probs))
    return np.mean(logits)


def get_mean_logit_mod_level(in_read):
    mod_mean_logit = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_probs = [this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]
        mod_mean_logit[this_mod] = get_mean_logit(this_mod_probs)
    return mod_mean_logit


base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'

# ds = 'WT'
# bam_file = os.path.join(base_dir, 'HEK293/WT_P2/chrALL.mAFiA.reads.bam')

# ds = 'M3KO'
# bam_file = os.path.join(base_dir, 'HEK293T_Mettl3_KO/merged/chrALL.mAFiA.reads.bam')

# ds = 'M3KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siMETTL3_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siTRUB1_input_merged/chrALL.mAFiA.reads.bam')

ds = 'TRUB1OE'
bam_file = os.path.join(base_dir, 'HEK293_TRUB1_OE/merged/chrALL.mAFiA.reads.bam')

img_out = '/home/adrian/img_out/single_read_cross_talk'
os.makedirs(img_out, exist_ok=True)

# gene_bed = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/gene.ensembl_havana.GRCh38.102.bed'
# df_gene = pd.read_csv(gene_bed, sep='\t')

single_read_mean_logit = []
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for this_read in tqdm(bam.fetch()):
        single_read_mean_logit.append(get_mean_logit_mod_level(this_read))

vec_m6A, vec_psi = np.vstack([
    (this_read_mean_logit['m6A'],  this_read_mean_logit['psi'])
    for this_read_mean_logit in single_read_mean_logit
    if ~np.isnan(this_read_mean_logit['m6A']) and ~np.isnan(this_read_mean_logit['psi'])
]).T

num_valid_reads = len(vec_m6A)

num_hh = np.sum((vec_m6A > 0) * (vec_psi > 0))
num_hl = np.sum((vec_m6A > 0) * (vec_psi < 0))
num_lh = np.sum((vec_m6A < 0) * (vec_psi > 0))
num_ll = np.sum((vec_m6A < 0) * (vec_psi < 0))
num_total = num_hh + num_hl + num_lh + num_ll

perc_hh = round(num_hh / num_total * 100, 2)
perc_hl = round(num_hl / num_total * 100, 2)
perc_lh = round(num_lh / num_total * 100, 2)
perc_ll = round(num_ll / num_total * 100, 2)

xy_max = 10
num_bins = 80
mat_z, edges_x, edges_y = np.histogram2d(vec_m6A, vec_psi, bins=num_bins, range=[[-xy_max, xy_max], [-xy_max, xy_max]])

plt.figure(figsize=(4, 4))
plt.axvline(x=0, c='red', ls='--')
plt.axhline(y=0, c='red', ls='--')
plt.imshow(np.log10(mat_z+1), extent=[-xy_max, xy_max, -xy_max, xy_max], origin='lower', vmax=2.5)
ax = plt.gca()
plt.text(0.1, 0.1, f'{perc_ll}%', c='r', fontsize=12, ha='left', va='bottom', transform=ax.transAxes)
plt.text(0.1, 0.9, f'{perc_lh}%', c='r', fontsize=12, ha='left', va='top', transform=ax.transAxes)
plt.text(0.9, 0.1, f'{perc_hl}%', c='r', fontsize=12, ha='right', va='bottom', transform=ax.transAxes)
plt.text(0.9, 0.9, f'{perc_hh}%', c='r', fontsize=12, ha='right', va='top', transform=ax.transAxes)
plt.xlim([-xy_max, xy_max])
plt.ylim([-xy_max, xy_max])
plt.xticks(np.linspace(-xy_max, xy_max, 5))
plt.yticks(np.linspace(-xy_max, xy_max, 5))
plt.xlabel(rf"$\langle$logit $p({{{dict_mod_display['m6A']}}})$$\rangle$ per read")
plt.ylabel(rf"$\langle$logit $p({{{dict_mod_display['psi']}}})$$\rangle$ per read")
plt.title(f'{ds}\n{num_valid_reads} valid reads')
plt.savefig(os.path.join(img_out, f'mean_logit_S_per_read_{ds}.png'), bbox_inches='tight')
plt.close('all')

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
