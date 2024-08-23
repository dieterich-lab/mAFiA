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

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

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


thresh_IntronDepth = 10
thresh_IRratio = 0.05
thresh_intronic_percentage = 0.1

ds = 'TAC'
cond = 'TAC_merged'
res_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/{cond}'
irf_file = os.path.join(res_dir, 'IRFinder/IRFinder-IR-nondir.txt')
bam_file = os.path.join(res_dir, 'chrALL.mAFiA.reads.bam')

img_out = '/home/adrian/img_out/IR_reads_avg_mod'
os.makedirs(img_out, exist_ok=True)

df_irf = pd.read_csv(irf_file, sep='\t', dtype={'Chr': str})
df_irf_thresh = df_irf[
    (df_irf['IntronDepth'] >= thresh_IntronDepth)
    * (df_irf['IRratio'] >= thresh_IRratio)
]
df_irf_clean = df_irf_thresh.loc[[this_name.split('/')[-1] == 'clean' for this_name in df_irf_thresh['Name']], :]

intronic_reads = []
exonic_reads = []
for _, this_row in tqdm(df_irf_clean.iterrows()):
    flag_required = 0 if this_row['Strand'] == '+' else 16
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for this_read in bam.fetch(this_row['Chr'], this_row['Start'], this_row['End']):
            if this_read.flag == flag_required:
                ref_to_read_pos = {ref_pos: read_pos for read_pos, ref_pos in this_read.get_aligned_pairs()}
                intronic_read_pos = [ref_to_read_pos.get(this_pos) for this_pos in range(this_row['Start'], this_row['End'])]
                intronic_percentage = len([x for x in intronic_read_pos if x is not None]) / len(intronic_read_pos)
                if intronic_percentage >= thresh_intronic_percentage:
                    intronic_reads.append(this_read)
                else:
                    exonic_reads.append(this_read)

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
plt.savefig(os.path.join(img_out, f'hist_{ds}_{cond}.png'), bbox_inches='tight')