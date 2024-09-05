import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import re
import numpy as np
from tqdm import tqdm
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
#######################################################################
import matplotlib.pyplot as plt


def get_n(in_df, in_chrom, in_chromStart, in_chromEnd, in_strand):
    sel_df = in_df[
        (in_df['chrom'] == in_chrom)
        * (in_df['chromStart'] >= in_chromStart)
        * (in_df['chromStart'] < in_chromEnd)
        * (in_df['strand'] == in_strand)
        ]
    return len(sel_df)


def get_mod_gene_log2fc(min_sites_per_gene=5):
    out_mod_gene_log2fc = {
        this_mod: {} for this_mod in mods
    }
    for this_mod in mods:
        print(this_mod)
        for this_gene in tqdm(df_gene['gene'].unique()):
            sub_df = df_gene[df_gene['gene'] == this_gene]
            this_chrom = sub_df['chrom'].unique()[0]
            this_chromStart = sub_df['chromStart'].min()
            this_chromEnd = sub_df['chromEnd'].max()
            this_strand = sub_df['strand'].unique()[0]
            cond_norm_mod_sites = []
            for this_cond in conditions:
                df_mod = cond_mod_thresh_dfs[this_cond][this_mod][thresholds[1]]
                n_mod = get_n(df_mod, this_chrom, this_chromStart, this_chromEnd, this_strand)
                df_covered = cond_mod_thresh_dfs[this_cond][this_mod][thresholds[0]]
                n_covered = get_n(df_covered, this_chrom, this_chromStart, this_chromEnd, this_strand)
                if n_mod >= min_sites_per_gene:
                    cond_norm_mod_sites.append(n_mod / n_covered)
                else:
                    cond_norm_mod_sites.append(None)
            if (cond_norm_mod_sites[0] is not None) and (cond_norm_mod_sites[1] is not None):
                this_log2fc = np.log2(cond_norm_mod_sites[1] / cond_norm_mod_sites[0])
                out_mod_gene_log2fc[this_mod][this_gene] = (cond_norm_mod_sites[0], cond_norm_mod_sites[1], this_log2fc)
                # print(this_gene, out_mod_gene_log2fc[this_mod][this_gene])
    return out_mod_gene_log2fc


bed6_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand'
]

gene_bed_fields = [
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

gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'
df_gene = pd.read_csv(gene_bed, sep='\t', names=gene_bed_fields)
df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR'
thresholds = ['0.0', '50.0']
mods = ['m6A', 'psi']

ds = 'TAC'
conditions = ['SHAM_merged', 'TAC_merged']

# ds = 'HFpEF'
# conditions = ['ctrl', 'HFpEF']

# ds = 'Diet'
# conditions = ['WT_CD_merged', 'WT_WD_merged']

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure3'
os.makedirs(img_out, exist_ok=True)

cond_mod_thresh_dfs = {
    this_cond: {
        this_mod: {} for this_mod in mods
    } for this_cond in conditions
}
for this_cond in conditions:
    for this_mod in mods:
        for this_thresh in thresholds:
            this_filename = os.path.join(res_dir, f'{ds}_{this_cond}_{this_mod}_modRatio{this_thresh}.bed')
            this_df = pd.read_csv(this_filename, sep='\t', names=bed6_fields)
            this_df['chrom'] = [this_chr.lstrip('chr') for this_chr in this_df['chrom']]
            cond_mod_thresh_dfs[this_cond][this_mod][this_thresh] = this_df


file_mod_gene_log2fc = os.path.join(img_out, f'df_mod_gene_log2fc_{ds}.tsv')

if os.path.exists(file_mod_gene_log2fc):
    df_mod_gene_log2fc = pd.read_csv(file_mod_gene_log2fc, sep='\t')
    mod_gene_log2fc = {this_mod: {} for this_mod in mods}
    for this_mod in mods:
        mod_df = df_mod_gene_log2fc[df_mod_gene_log2fc['mod']==this_mod]
        for _, this_row in mod_df.iterrows():
            mod_gene_log2fc[this_mod][this_row['gene']] = (
                this_row['mod_ratio_ctrl'], this_row['mod_ratio_HFpEF'], this_row['log2fc']
            )
else:
    mod_gene_log2fc = get_mod_gene_log2fc()
    df_mod_gene_log2fc = pd.DataFrame(
        [[k, kk, *vv] for k, v in mod_gene_log2fc.items() for kk, vv in v.items()],
        columns=['mod', 'gene', f'mod_ratio_{conditions[0]}', f'mod_ratio_{conditions[1]}', 'log2fc']
    )
    df_mod_gene_log2fc.to_csv(file_mod_gene_log2fc, sep='\t', index=False, float_format='%.3f')


max_display_genes = 4

plt.figure(figsize=(5*cm, 2*cm))
for mod_ind, this_mod in enumerate(mods):
    top_log2fc_genes = [item[0]
                        for item in sorted(mod_gene_log2fc[this_mod].items(),
                                           key=lambda item: item[1][2], reverse=True)[:max_display_genes]
                        if item[1][2] > 0.1][::-1]
    bottom_log2fc_genes = [item[0]
                           for item in sorted(mod_gene_log2fc[this_mod].items(),
                                              key=lambda item: item[1][2])[:max_display_genes]
                           if item[1][2] < -0.1]
    sorted_genes = bottom_log2fc_genes + top_log2fc_genes

    plt.subplot(1, 2, mod_ind+1)
    if mod_ind == 1:
        plt.gca().yaxis.tick_right()
    plt.barh(range(len(sorted_genes)), [mod_gene_log2fc[this_mod][this_gene][2] for this_gene in sorted_genes])
    plt.xlim([-1, 1])
    plt.xticks(np.arange(-1, 1.01, 0.5))
    plt.yticks(range(len(sorted_genes)), sorted_genes)
    plt.axvline(x=0, c='gray', ls='--', lw=0.5)
plt.savefig(os.path.join(img_out, f'log2fc_num_mod_sites_{ds}.{FMT}'), **fig_kwargs)
    