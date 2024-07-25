import os
import pandas as pd
import numpy as np
import pysam
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import kstest

day = 'day21'
ks_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/KS_test/annotated_sites_{day}.tsv'
polyA_sham_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA/SHAM_{day}_read2polyA_length.txt'
polyA_tac_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA/TAC_{day}_read2polyA_length.txt'
sham_bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_{day}/chrALL.mAFiA.reads.bam'
tac_bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_{day}/chrALL.mAFiA.reads.bam'
gene_bed_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'

img_out = '/home/adrian/img_out/polyA_len'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'feature_type',
    'other',
    'description'
]

def get_gene_read_ids(in_bam_file, in_chrom, in_strand, in_featureStart, in_featureEnd):
    read_ids = []
    flag_require = 0 if in_strand == '+' else 16
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in bam.fetch(in_chrom, in_featureStart, in_featureEnd):
            if read.flag == flag_require:
                read_ids.append(read.qname)
    return read_ids


def collect_polyA_len(in_df, in_read_ids, thresh_len=0):
    out_len = []
    for this_read_id in in_read_ids:
        sel_df = in_df[in_df['read_id']==this_read_id]
        if len(sel_df):
            val = sel_df['polyA_len'].values[0]
            if val >= thresh_len:
                out_len.append(val)
    return out_len


df_ks = pd.read_csv(ks_file, sep='\t')
df_ks_sel = df_ks[
    (df_ks['feature']=='three_prime_utr')
    * (np.abs(df_ks['delta']) >= 10.0)
    ]
target_genes = df_ks_sel['gene_name'].unique()
df_polyA_sham = pd.read_csv(polyA_sham_file, sep='\t', names=['read_id', 'polyA_len'])
df_polyA_tac = pd.read_csv(polyA_tac_file, sep='\t', names=['read_id', 'polyA_len'])

plt.figure(figsize=(5, 4))
plt.hist(df_polyA_sham['polyA_len'], range=[0, 400], bins=50, alpha=0.5, label='SHAM')
plt.hist(df_polyA_tac['polyA_len'], range=[0, 400], bins=50, alpha=0.5, label='TAC')
plt.legend(loc='upper right')
plt.xlabel('polyA length', fontsize=12)
plt.ylabel('Read count', fontsize=12)
plt.title(day, fontsize=15)
plt.savefig(os.path.join(img_out, 'hist_polyA_length_overall.png'), bbox_inches='tight')

gene_bed = pd.read_csv(gene_bed_file, sep='\t', names=bed_fields)
gene_bed['gene_name'] = [field.split(' ')[1].strip('\"')
                         for desc in gene_bed['description']
                         for field in desc.split('; ') if field.split(' ')[0] == 'gene_name']

delta_polyA_lens = {}
for this_gene in tqdm(target_genes):
    df_gene = df_ks_sel[df_ks_sel['gene_name'] == this_gene]
    chrom, strand = df_gene[['chrom', 'strand']].values[0]
    geneStart, geneEnd = gene_bed[gene_bed['gene_name']==this_gene][['chromStart', 'chromEnd']].values[0]
    sham_read_ids = get_gene_read_ids(sham_bam_file, chrom, strand, geneStart, geneEnd)
    tac_read_ids = get_gene_read_ids(tac_bam_file, chrom, strand, geneStart, geneEnd)

    sham_polyA_len = collect_polyA_len(df_polyA_sham, sham_read_ids)
    tac_polyA_len = collect_polyA_len(df_polyA_tac, tac_read_ids)
    ks_stat, pval = kstest(sham_polyA_len, tac_polyA_len)

    delta_polyA_lens[this_gene] = (list(zip(df_gene['name'].values, df_gene['delta'].values)),
                                   np.mean(sham_polyA_len), np.mean(tac_polyA_len),
                                   len(sham_polyA_len), len(tac_polyA_len), ks_stat, pval)

delta_mod_polyA_len = {}
for k, v in delta_polyA_lens.items():
    if v[6]>=0.3:
        continue
    polyA_delta = v[2] - v[1]
    mod_mean_deltas = {}
    for this_mod in mods:
        deltas = [delta for mod, delta in v[0] if mod==this_mod]
        if len(deltas):
            mod_mean_deltas[this_mod] = np.mean(deltas)
            # mod_mean_deltas[this_mod] = deltas[np.argmax(np.abs(deltas))]
    if len(mod_mean_deltas.values()):
        delta_mod_polyA_len[k] = (mod_mean_deltas, polyA_delta)

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    # vec_x, vec_y = np.vstack([v[this_mod] for v in delta_mod_polyA_len.values()]).T
    vec_x, vec_y = np.vstack([(v[0].get(this_mod, np.nan), v[1]) for v in delta_mod_polyA_len.values()]).T
    plt.scatter(vec_x, vec_y)
    plt.xlim([-30, 30])
    plt.ylim([-100, 100])
    plt.xlabel(rf'$\Delta S_{{{dict_mod_display[this_mod]}}}$')
    plt.ylabel(r'$\Delta L(polyA)$ (bps)')
plt.suptitle(f'{day}', fontsize=15)

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    # vec_x, vec_y = np.vstack([v[this_mod] for v in delta_mod_polyA_len.values()]).T
    vec_x, vec_y = np.vstack([(v[0].get(this_mod, np.nan), v[1]) for v in delta_mod_polyA_len.values()]).T
    down_changes = vec_y[vec_x<0]
    up_changes = vec_y[vec_x>0]
    plt.violinplot([down_changes, up_changes], range(2), showmeans=True)
    # down_change = np.mean(vec_y[vec_x<0])
    # up_change = np.mean(vec_y[vec_x>0])
    # plt.bar(range(2), [down_change, up_change], width=0.4)
    plt.xticks(range(2), ['down', 'up'])
    plt.xlim([-0.5, 1.5])
    # plt.ylim([-2.5, 2.5])
    # plt.scatter(vec_x, vec_y)
    # plt.xlim([-30, 30])
    # plt.ylim([-100, 100])
    plt.axhline(y=0, c='gray')
    plt.xlabel(rf'$\Delta S_{{{dict_mod_display[this_mod]}}}$')
    plt.ylabel(r'$\Delta L(polyA)$ (bps)')
plt.suptitle(f'{day}', fontsize=15)

# plot_genes = [
#     'Ube2d3',
#     'Pdia3',
#     'Csnk1a1',
#     'H2-K1'
# ]
# for this_gene in plot_genes:
#     df_gene = df_ks_sel[df_ks_sel['gene_name'] == this_gene]
#     chrom, strand, featureStart, featureEnd = df_gene[['chrom', 'strand', 'featureStart', 'featureEnd']].values[0]
#     sham_read_ids = get_gene_read_ids(sham_bam_file, chrom, strand, featureStart, featureEnd)
#     tac_read_ids = get_gene_read_ids(tac_bam_file, chrom, strand, featureStart, featureEnd)
#     sham_polyA_len = collect_polyA_len(df_polyA_sham, sham_read_ids)
#     tac_polyA_len = collect_polyA_len(df_polyA_tac, tac_read_ids)
#     ks_stat, pval = kstest(sham_polyA_len, tac_polyA_len)
#
#     plt.figure(figsize=(5, 5))
#     plt.hist(sham_polyA_len, bins=50, range=[0, 400], fc='b', alpha=0.5, label='SHAM')
#     plt.hist(tac_polyA_len, bins=50, range=[0, 400], fc='r', alpha=0.5, label='TAC')
#     plt.legend(loc='upper left')
#     plt.xlabel('polyA length (bps)', fontsize=12)
#     plt.ylabel('Counts', fontsize=12)
#     plt.title(f'{this_gene}, {day}', fontsize=15)
