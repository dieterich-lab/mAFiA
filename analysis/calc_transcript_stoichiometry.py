import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import pysam
import numpy as np
import re
from scipy.stats import kstest
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'
out_dir = os.path.join(res_dir, 'transcript_stoichiometry')
os.makedirs(out_dir, exist_ok=True)

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']

gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'frame',
    'description'
]
df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
# df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}


def get_mod_probs_per_gene(in_bam, in_df):
    gene_mod_probs = {}
    for _, this_row in tqdm(in_df.iterrows()):
        chrom, chromStart, chromEnd, strand, gene = this_row[['chrom', 'chromStart', 'chromEnd', 'strand', 'gene']]
        gene_mod_probs[gene] = {mod: [] for mod in mods}
        flag_required = 0 if strand == '+' else 16
        for read in in_bam.fetch(chrom, chromStart, chromEnd):
            if read.flag == flag_required:
                for mod in mods:
                    gene_mod_probs[gene][mod].extend([pair[1] for pair in read.modified_bases.get(('N', 0, dict_mod_code[mod]), [])])
    return gene_mod_probs


def get_num_nts_log2fc_pval(in_mod_probs_0, in_mod_probs_1):
    num_nts_0 = len(in_mod_probs_0)
    num_nts_1 = len(in_mod_probs_1)
    if (num_nts_0 > 0) and (num_nts_1 > 0):
        stoichio_0 = np.mean(np.array(in_mod_probs_0) >= 128)
        stoichio_1 = np.mean(np.array(in_mod_probs_1) >= 128)
        log2fc_stoichio = np.log2(stoichio_1 / stoichio_0)
        ks_stat, pval = kstest(in_mod_probs_0, in_mod_probs_1)
        return num_nts_0, num_nts_1, stoichio_0, stoichio_1, log2fc_stoichio, pval
    else:
        return None

with pysam.AlignmentFile(os.path.join(res_dir, ds, conditions[0], 'chrALL.mAFiA.reads.bam'), 'rb') as bam0:
    mod_probs_0 = get_mod_probs_per_gene(bam0, df_gene)

with pysam.AlignmentFile(os.path.join(res_dir, ds, conditions[1], 'chrALL.mAFiA.reads.bam'), 'rb') as bam1:
    mod_probs_1 = get_mod_probs_per_gene(bam1, df_gene)


gene_mod_log2fc = {}
for this_gene in tqdm(df_gene['gene']):
    gene_mod_log2fc[this_gene] = {}
    for this_mod in mods:
        gene_mod_log2fc[this_gene][this_mod] = get_num_nts_log2fc_pval(mod_probs_0[this_gene][this_mod],
                                                                       mod_probs_1[this_gene][this_mod])
gene_mod_log2fc = {k: v for k, v in gene_mod_log2fc.items() if any([vv is not None for vv in v.values()])}

df_out = pd.DataFrame([[k] + [kk] + list(vv) for k, v in gene_mod_log2fc.items() for kk, vv in v.items() if vv is not None],
                      columns=['gene',
                               'mod',
                               f'num_nts_{conditions[0]}',
                               f'num_nts_{conditions[1]}',
                               f'stoichiometry_{conditions[0]}',
                               f'stoichiometry_{conditions[1]}',
                               'log2fc',
                               'pval'
                               ]
                      )

df_out.to_csv(os.path.join(out_dir, f'{ds}_{conditions[1]}_vs_{conditions[0]}.tsv'),
              sep='\t', index=False, float_format='%.3f')

### volcano plot ###
xmax = 4
ymax = 10
thresh_log_pval = 4
thresh_log2fc = 0.1
reg_genes = {mod: {} for mod in mods}
plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    this_mod_df = df_out[df_out['mod'] == this_mod]
    vec_gene, vec_log2fc, vec_pval = this_mod_df[['gene', 'log2fc', 'pval']].values.T
    vec_log_pval = -np.log10(np.float64(vec_pval))

    mask_up = (vec_log_pval >= thresh_log_pval) * (vec_log2fc > thresh_log2fc)
    mask_down = (vec_log_pval >= thresh_log_pval) * (vec_log2fc < -thresh_log2fc)

    plt.subplot(1, 2, mod_ind+1)
    plt.scatter(vec_log2fc, vec_log_pval, s=1, c='gray', alpha=0.5)
    plt.scatter(vec_log2fc[mask_up], vec_log_pval[mask_up], s=3, c='red')
    plt.scatter(vec_log2fc[mask_down], vec_log_pval[mask_down], s=3, c='blue')
    plt.xlim([-xmax, xmax])
    plt.ylim([0, ymax])
    plt.xlabel(f'log2fc transcript ${{{dict_mod_display[this_mod]}}}$')
    plt.ylabel('-log10(pval)')

    reg_genes[this_mod]['up'] = vec_gene[mask_up]
    reg_genes[this_mod]['down'] = vec_gene[mask_down]

    # reg_genes[this_mod]['up'] = [this_gene for this_gene in vec_gene[mask_up] if this_gene[:2]!='Gm']
    # reg_genes[this_mod]['down'] = [this_gene for this_gene in vec_gene[mask_down] if this_gene[:2]!='Gm']

    # plt.text(-xmax+1, 1, '\n'.join(reg_genes[this_mod]['down']), c='blue')
    # plt.text(xmax-2, 1, '\n'.join(reg_genes[this_mod]['up']), c='red')

plt.savefig(os.path.join(out_dir, f'volcano_plot_logPval{thresh_log_pval}_log2fc{thresh_log2fc}.png'),
            bbox_inches='tight')

with open(os.path.join(out_dir, f'reg_genes_logPval{thresh_log_pval}_log2fc{thresh_log2fc}.txt'), 'w') as f_out:
    for this_mod in mods:
        for this_reg in ['up', 'down']:
            f_out.write(f'{this_mod}_{this_reg}\n')
            f_out.write('\n'.join(reg_genes[this_mod][this_reg]))
            f_out.write('\n\n')