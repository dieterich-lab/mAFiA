import pandas as pd
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import kstest
import numpy as np
from scipy.stats import trim_mean, ttest_ind

this_day = 'day56'

thresh_confidence = 0.0
thresh_coverage = 20
thresh_pval = 1.0
thresh_num_sites = 1
mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
mod_color = {
    'm6A': 'r',
    'psi': 'b'
}
condition_color = {
    'SHAM': 'blue',
    'TAC': 'red'
}

img_out = '/home/adrian/img_out/log2fc_per_gene'
os.makedirs(img_out, exist_ok=True)

gene_annotation = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.gtf'
df_gene_gtf = pd.read_csv(gene_annotation, sep='\t',
            names=[
                'seqname',
                'source',
                'feature',
                'start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute',
            ]
            )
df_gene_bed = []
for _, this_row in df_gene_gtf.iterrows():
    chrom = this_row['seqname']
    chromStart = this_row['start'] - 1
    chromEnd = this_row['end'] - 1
    gene_name = [field.split(' ')[1] for field in this_row['attribute'].split('; ') if field.split(' ')[0]=='gene_name']
    name = gene_name[0].strip('\"')
    score = this_row['score']
    strand = this_row['strand']
    df_gene_bed.append([chrom, chromStart, chromEnd, name, score, strand])
df_gene_bed = pd.DataFrame(df_gene_bed, columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

# proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis.xlsx'
# df_proteome = pd.read_excel(proteome_file)
proteome_file = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis_benjamini_hochberg.tsv'
df_proteome = pd.read_csv(proteome_file, sep='\t')
df_proteome = df_proteome[
    (df_proteome['main.comparison']=='TAC vs Sham')
    * (df_proteome['tp']==('Day ' + this_day.lstrip('day')))
    * (df_proteome['p.adj']<=thresh_pval)
    ]
gene_protein_log2fc = {}
for this_gene in df_proteome['SYMBOL'].unique():
    sub_df = df_proteome[df_proteome['SYMBOL']==this_gene]
    gene_protein_log2fc[this_gene] = (sub_df['logFC'].values[0], sub_df['p.adj'].values[0])

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'
conditions = ['SHAM', 'TAC']
condition_days = [f'{this_cond}_{this_day}' for this_cond in conditions]

log2fc_file = os.path.join(img_out, f'gene_log2fc_{this_day}.tsv')

### merge SHAM and TAC ###
dfs = {}
for this_cond in conditions:
    df_in = pd.read_csv(os.path.join(res_dir, f'{this_cond}_{this_day}', 'chrALL.mAFiA.sites.bed'), sep='\t', dtype={'chrom': 'str'})
    dfs[f'{this_cond}_{this_day}'] = df_in[
        (df_in['confidence']>=thresh_confidence)
        * (df_in['coverage']>=thresh_coverage)
    ]
df_merged = pd.merge(*list(dfs.values()),
                     on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                     suffixes=[f'_{this_key}' for this_key in dfs.keys()])


def collect_gene_mod_values(gene_list):
    gene_mod_values = {}
    for _gene in tqdm(gene_list):
        gene_subdf = df_gene_bed[df_gene_bed['name'] == _gene]
        _chrom = gene_subdf['chrom'].unique()[0]
        _chromStart = gene_subdf['chromStart'].min()
        _chromEnd = gene_subdf['chromEnd'].max()
        _strand = gene_subdf['strand'].unique()[0]
        df_gene_mod_sites = df_merged[
            (df_merged['chrom'] == _chrom)
            * (df_merged['chromStart'] >= _chromStart)
            * (df_merged['chromEnd']<_chromEnd)
            * (df_merged['strand']==_strand)
        ]
        this_gene_mod_values = {
            my_mod: {
                this_cond_day: df_gene_mod_sites[df_gene_mod_sites['name'] == my_mod][f'modRatio_{this_cond_day}'].values
                for this_cond_day in condition_days
            }
            for my_mod in ['m6A', 'psi']
        }

        if np.concatenate([this_cond_vals for this_mod_vals in this_gene_mod_values.values() for this_cond_vals in this_mod_vals.values()]).any():
            gene_mod_values[_gene] = this_gene_mod_values
    return gene_mod_values


def calc_log2fc_stats():
    gene_mod_values = collect_gene_mod_values(df_gene_bed['name'].unique())

    ### log2FC ###
    metric = 'mean'
    my_gene_log2fc_stats = {}
    for gene, mod_ratios in gene_mod_values.items():
            my_gene_log2fc_stats[gene] = {}
            for _mod in ['m6A', 'psi']:
                if np.concatenate(list(mod_ratios[_mod].values())).any():
                    num_sites = len(mod_ratios[_mod][f'SHAM_{this_day}'])
                    sham_vals = mod_ratios[_mod][f'SHAM_{this_day}']
                    tac_vals = mod_ratios[_mod][f'TAC_{this_day}']
                    if metric=='mean':
                        sham_metric = np.mean(sham_vals)
                        tac_metric = np.mean(tac_vals)
                    elif metric=='trimmed_mean':
                        sham_metric = trim_mean(sham_vals, 0.1)
                        tac_metric = trim_mean(tac_vals, 0.1)
                    elif metric=='median':
                        sham_metric = np.median(sham_vals)
                        tac_metric = np.median(tac_vals)
                    elif metric=='75percentile':
                        sham_metric = np.percentile(sham_vals, q=75)
                        tac_metric = np.percentile(tac_vals, q=75)
                    elif metric=='90percentile':
                        sham_metric = np.percentile(sham_vals, q=90)
                        tac_metric = np.percentile(tac_vals, q=90)
                    elif metric=='max':
                        sham_metric = np.max(sham_vals)
                        tac_metric = np.max(tac_vals)
                    tac_metric += 0.000001
                    sham_metric += 0.000001
                    log2fc = np.log2(tac_metric / sham_metric)
                    # t_stat, p_val = ttest_ind(sham_vals, tac_vals, equal_var=False)
                    ks_stat, p_val = kstest(sham_vals, tac_vals)
                    my_gene_log2fc_stats[gene][_mod] = (num_sites, sham_metric, tac_metric, log2fc, p_val)

    return my_gene_log2fc_stats

if os.path.exists(log2fc_file):
    df_log2fc = pd.read_csv(log2fc_file, sep='\t')
    gene_log2fc_stats = {}
    for _, this_row in df_log2fc.iterrows():
        if this_row['gene'] not in gene_log2fc_stats.keys():
            gene_log2fc_stats[this_row['gene']] = {}
        gene_log2fc_stats[this_row['gene']][this_row['mod']] = (this_row['num_mod_sites'],
                                                                this_row['SHAM'],
                                                                this_row['TAC'],
                                                                this_row['log2fc_mod'],
                                                                this_row['p_val_mod']
                                                                )
else:
    gene_log2fc_stats = calc_log2fc_stats()

### distribution of num. sites per gene ###
plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    plt.hist([stats.get(this_mod)[0] for gene, stats in gene_log2fc_stats.items() if stats.get(this_mod)], bins=50, range=[0, 100])
    plt.xlabel('Num. comparable sites in gene', fontsize=12)
    plt.ylabel('Gene count', fontsize=12)
plt.suptitle(f'{this_day}', fontsize=15)
plt.savefig(os.path.join(img_out, 'hist_num_mod_sites_per_gene.png'), bbox_inches='tight')

### all ###
max_range = 1.5
gene_log2fc_m6A = {
    this_gene: this_gene_log2fc.get('m6A')[3]
    for this_gene, this_gene_log2fc in gene_log2fc_stats.items()
    if (this_gene_log2fc.get('m6A') is not None)
       and (this_gene_log2fc.get('m6A')[4] <= thresh_pval)
       and (this_gene_log2fc.get('m6A')[0] >= thresh_num_sites)
}
gene_log2fc_psi = {
    this_gene: this_gene_log2fc.get('psi')[3]
    for this_gene, this_gene_log2fc in gene_log2fc_stats.items()
    if (this_gene_log2fc.get('psi') is not None)
       and (this_gene_log2fc.get('psi')[4] <= thresh_pval)
       and (this_gene_log2fc.get('psi')[0] >= thresh_num_sites)
}
all_log2fc = {
    'm6A': gene_log2fc_m6A.values(),
    'psi': gene_log2fc_psi.values()
}
# plt.figure(figsize=(10, 5))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     plt.hist(all_log2fc[this_mod], bins=50, range=[-max_range, max_range], alpha=0.5)
#     plt.axvline(x=0, c='g')
#     plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
#     plt.ylabel('Gene counts', fontsize=12)
#     plt.title(f'N={len(all_log2fc[this_mod])}', fontsize=15)


### output table ###
df_log2fc_stats = []
for this_gene, this_log2fc_stats in gene_log2fc_stats.items():
    this_protein_log2fc_stats = gene_protein_log2fc.get(this_gene)
    for this_mod in mods:
        this_mod_log2fc_stats = this_log2fc_stats.get(this_mod)
        if (this_mod_log2fc_stats is not None) and (this_protein_log2fc_stats is not None):
            df_log2fc_stats.append([this_gene, this_mod] + list(this_mod_log2fc_stats) + list(this_protein_log2fc_stats))
if not os.path.exists(log2fc_file):
    df_log2fc_stats = pd.DataFrame(df_log2fc_stats,
                                   columns=[
                                       'gene',
                                       'mod',
                                       'num_mod_sites',
                                       'SHAM',
                                       'TAC',
                                       'log2fc_mod',
                                       'p_val_mod',
                                       'log2fc_protein',
                                       'p_val_protein'
                                   ])
    df_log2fc_stats.to_csv(log2fc_file, sep='\t', index=False, float_format='%.3f')

### common gene and thresholded by min. sites ###
common_genes, common_log2fc_m6A, common_log2fc_psi = np.vstack(
    [(gene, val['m6A'][3], val['psi'][3]) for gene, val in gene_log2fc_stats.items()
     if (set(val.keys()).issuperset((set(['m6A', 'psi']))))
     and (val['m6A'][4] <= thresh_pval) and (val['psi'][4] <= thresh_pval)
     and (gene in gene_protein_log2fc.keys())
     ]
).T
common_log2fc_protein = [gene_protein_log2fc.get(this_gene)[0] if this_gene in gene_protein_log2fc.keys() else np.nan for this_gene in common_genes]
common_log2fc = {
    'm6A': np.float64(common_log2fc_m6A),
    'psi': np.float64(common_log2fc_psi),
    'protein': np.float64(common_log2fc_protein)
}

max_log2fc = 2.0
sep_log2fc = 1.0
log2fc_ranges = np.round(np.arange(-max_log2fc, max_log2fc+sep_log2fc, sep_log2fc), 1)

# colors = matplotlib.cm.tab20(range(len(log2fc_ranges)-1))
# plt.figure(figsize=(5, 5))
# for i in range(len(log2fc_ranges)-1):
#     range_start = log2fc_ranges[i]
#     range_end = log2fc_ranges[i+1]
#     mask = (common_log2fc['protein']>=range_start) * (common_log2fc['protein']<range_end)
#     plt.scatter(common_log2fc['m6A'][mask], common_log2fc['psi'][mask], s=1, color=colors[i], label=f'[{range_start}, {range_end})')
# plt.xlim([-max_log2fc, max_log2fc])
# plt.ylim([-max_log2fc, max_log2fc])
# plt.legend()

log2fc_combos = [
    r'$\downarrow \downarrow$',
    r'$\downarrow \uparrow$',
    r'$\uparrow \downarrow$',
    r'$\uparrow \uparrow$',
]

all_binned_protein_log2fc = []
all_binned_protein_log2fc.append(common_log2fc['protein'][(common_log2fc['m6A'] < 0) * (common_log2fc['psi'] < 0)])
all_binned_protein_log2fc.append(common_log2fc['protein'][(common_log2fc['m6A'] < 0) * (common_log2fc['psi'] >= 0)])
all_binned_protein_log2fc.append(common_log2fc['protein'][(common_log2fc['m6A'] >= 0) * (common_log2fc['psi'] < 0)])
all_binned_protein_log2fc.append(common_log2fc['protein'][(common_log2fc['m6A'] >= 0) * (common_log2fc['psi'] >= 0)])

mean_protein_log2fc = [np.nanmean(this_bin) for this_bin in all_binned_protein_log2fc]

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.violinplot(all_binned_protein_log2fc, range(len(log2fc_combos)))
plt.xticks(range(len(log2fc_combos)), log2fc_combos)
plt.ylim([-1.5, 1.5])
plt.xlabel(f"$log_{2}FC ({{{dict_mod_display['m6A']}}}, {{{dict_mod_display['psi']}}})$", fontsize=12)
plt.ylabel('$log_{2}FC (protein)$', fontsize=12)
plt.subplot(1, 2, 2)
plt.bar(range(len(log2fc_combos)), mean_protein_log2fc)
plt.xticks(range(len(log2fc_combos)), log2fc_combos)
plt.ylim([-0.1, 0.1])
plt.xlabel(f"$log_{2}FC ({{{dict_mod_display['m6A']}}}, {{{dict_mod_display['psi']}}})$", fontsize=12)
plt.ylabel(r'$\langle log_{2}FC (protein) \rangle$', fontsize=12)
plt.suptitle(this_day, fontsize=15)
plt.savefig(os.path.join(img_out, f'log2fc_protein_by_mod_changes_{this_day}.png'), bbox_inches='tight')

# all_binned_protein_log2fc = [[[] for j in range(len(log2fc_ranges)-1)] for i in range(len(log2fc_ranges)-1)]
# for row_i in range(len(log2fc_ranges)-1):
#     row_start = log2fc_ranges[row_i]
#     row_end = log2fc_ranges[row_i+1]
#     for col_i in range(len(log2fc_ranges)-1):
#         col_start = log2fc_ranges[col_i]
#         col_end = log2fc_ranges[col_i + 1]
#         mask = (common_log2fc['m6A']>=col_start) * (common_log2fc['m6A']<col_end) \
#                * (common_log2fc['psi']>=row_start) * (common_log2fc['psi']<row_end)
#         all_binned_protein_log2fc[row_i][col_i] = common_log2fc['protein'][mask]
# mean_binned_protein_log2fc = np.nan_to_num([[np.nanmean(this_col) for this_col in this_row] for this_row in all_binned_protein_log2fc])
#
# vrange = 0.2
#
# plt.figure(figsize=(5, 5))
# im = plt.imshow(mean_binned_protein_log2fc, cmap='cool', origin='lower', vmin=-vrange, vmax=vrange)
# plt.xticks(np.arange(len(log2fc_ranges))-0.5, log2fc_ranges)
# plt.yticks(np.arange(len(log2fc_ranges))-0.5, log2fc_ranges)
# plt.xlabel(f"$log_{2}FC ({{{dict_mod_display['m6A']}}})$", fontsize=12)
# plt.ylabel(f"$log_{2}FC ({{{dict_mod_display['psi']}}})$", fontsize=12)
# cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
# cbar.set_ticks(np.linspace(-vrange, vrange, 5))
# plt.title(f'{this_day}\n$log_{2}FC$ (protein)', fontsize=15)
# plt.savefig(os.path.join(img_out, f'hist2d_log2fc_protein_vs_mods_{this_day}.png'), bbox_inches='tight')

# plt.figure(figsize=(10, 5))
# for mod_ind, this_mod in enumerate(mods):
#     that_mod = [m for m in mods if m!=this_mod][0]
#     plt.subplot(1, 2, mod_ind+1)
#     plt.hist(common_log2fc[this_mod][common_log2fc[that_mod]>=0], bins=50, range=[-max_range, max_range], fc='b', label=f'${dict_mod_display[that_mod]}$ up', alpha=0.5)
#     plt.hist(common_log2fc[this_mod][common_log2fc[that_mod]<0], bins=50, range=[-max_range, max_range], fc='r', label=f'${dict_mod_display[that_mod]}$ down', alpha=0.5)
#     plt.axvline(x=0, c='g')
#     plt.legend(loc='upper left')
#     plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
#     plt.ylabel('Gene counts', fontsize=12)
#     plt.title(f'N={len(common_log2fc[this_mod])}', fontsize=15)


### set up bins ###
max_bin = max_log2fc
# bin_edges = np.linspace(-max_bin, max_bin, 2*max_bin+1)
bin_edges = log2fc_ranges
bin_width = bin_edges[1] - bin_edges[0]
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

### m6A vs psi ###
plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    that_mod = [mod for mod in mods if mod!=this_mod][0]
    plt.subplot(1, 2, mod_ind+1)
    binned_that_mod_log2fc = []
    bin_counts = []
    for i in range(len(bin_edges)-1):
        bin_start = bin_edges[i]
        bin_end = bin_edges[i+1]
        binned_that_mod = common_log2fc[that_mod][
            (common_log2fc[this_mod] >= bin_start)
            * (common_log2fc[this_mod] < bin_end)
            ]
        bin_counts.append(len(binned_that_mod))
        binned_that_mod_log2fc.append(np.mean(binned_that_mod))
    total_counts = np.sum(bin_counts)
    plt.axhline(y=0, c='gray')
    plt.bar(bin_centers, binned_that_mod_log2fc, width=bin_width*0.8)
    for this_bin_center, this_bin_count, that_mod_log2fc in zip(bin_centers, bin_counts, binned_that_mod_log2fc):
        if not np.isnan(that_mod_log2fc):
            plt.text(this_bin_center, that_mod_log2fc + np.sign(that_mod_log2fc) * 0.1, this_bin_count,
                 horizontalalignment='center', verticalalignment='center')
    # plt.title(f'N={total_counts}')
    plt.xlim([-max_bin, max_bin])
    max_y = np.ceil(np.nanmax(np.abs(binned_that_mod_log2fc)))
    plt.ylim([-max_y, max_y])
    plt.xticks(bin_edges)
    plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
    plt.ylabel(f'$log_{2}FC ({dict_mod_display[that_mod]})$', fontsize=12)
plt.suptitle(f'{this_day}', fontsize=15)
plt.savefig(os.path.join(img_out, f'log2fc_m6A_vs_psi_{this_day}.png'), bbox_inches='tight')


### protein vs mods ###
common_genes_m6A = list(set(gene_log2fc_m6A.keys()).intersection(set(gene_protein_log2fc.keys())))
vec_log2fc_m6A, vec_log2fc_m6A_protein = np.vstack(
    [(gene_log2fc_m6A[this_gene], gene_protein_log2fc[this_gene][0]) for this_gene in common_genes_m6A]
).T
common_genes_psi = list(set(gene_log2fc_psi.keys()).intersection(set(gene_protein_log2fc.keys())))
vec_log2fc_psi, vec_log2fc_psi_protein = np.vstack(
    [(gene_log2fc_psi[this_gene], gene_protein_log2fc[this_gene][0]) for this_gene in common_genes_psi]
).T
vec_log2fc_mod = {}
vec_log2fc_mod['m6A'] = vec_log2fc_m6A
vec_log2fc_mod['psi'] = vec_log2fc_psi
vec_log2fc_protein = {}
vec_log2fc_protein['m6A'] = vec_log2fc_m6A_protein
vec_log2fc_protein['psi'] = vec_log2fc_psi_protein

# plt.figure(figsize=(10, 4))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     plt.scatter(vec_log2fc_mod[this_mod], vec_log2fc_protein[this_mod], s=2)
#     plt.xlim([-3, 3])
#     plt.ylim([-2, 2])
#     plt.axvline(x=0, c='r', ls='--')
#     plt.axhline(y=0, c='r', ls='--')
#     plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
#     plt.ylabel('$log_{2}FC (protein)$', fontsize=12)
# plt.suptitle(f'{this_day}', fontsize=15)

def plot_gene_mod_values(this_gene, this_gene_mod_vals, this_gene_log2fc_stats):
    plt.figure(figsize=(10, 5))
    for _ind, _mod in enumerate(mods):
        plt.subplot(1, 2, _ind+1)
        for cond_ind, this_cond in enumerate(conditions):
            violin_parts = plt.violinplot(this_gene_mod_vals[_mod][f'{this_cond}_{this_day}'],
                                          [cond_ind], showmeans=True, widths=0.5)
            for pc in violin_parts['bodies']:
                pc.set_facecolor(condition_color[this_cond])
                pc.set_edgecolor(condition_color[this_cond])
            violin_parts['cmaxes'].set_edgecolor(condition_color[this_cond])
            violin_parts['cmeans'].set_edgecolor(condition_color[this_cond])
            violin_parts['cmins'].set_edgecolor(condition_color[this_cond])
            violin_parts['cbars'].set_edgecolor(condition_color[this_cond])
        log2fc_stats = this_gene_log2fc_stats[_mod]
        plt.title(f'N={log2fc_stats[0]}, log2fc={log2fc_stats[3]}, pval={log2fc_stats[4]}')
        plt.xticks(range(len(conditions)), conditions)
        plt.xlim([-0.5, 1.5])
        plt.ylim([0, 100])
        plt.ylabel(f'$S_{{{dict_mod_display[_mod]}}}$', fontsize=12)
    plt.suptitle(f'{this_gene}\n{this_day}', fontsize=15)
    plt.savefig(os.path.join(img_out, f'{this_gene}_{this_day}.png'), bbox_inches='tight')

# sel_genes = np.array(common_genes_m6A)[(vec_log2fc_m6A>=1)*(vec_log2fc_m6A<2)]
# # sel_genes = ['Coa5']
# sel_gene_mod_vals = collect_gene_mod_values(sel_genes)
# for gene, mod_vals in sel_gene_mod_vals.items():
#     plot_gene_mod_values(gene, mod_vals, gene_log2fc_stats[gene])

### aggregate protein changes binned by mod changes ###
max_bin = 2
bin_edges = np.linspace(-max_bin, max_bin, 5)
bin_width = bin_edges[1] - bin_edges[0]
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
max_y = 1.5

mod_means = {}
plt.figure(figsize=(15, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 3, mod_ind+1)
    binned_protein_log2fc = []
    bin_counts = []
    for i in range(len(bin_edges)-1):
        bin_start = bin_edges[i]
        bin_end = bin_edges[i+1]
        binned_protein = vec_log2fc_protein[this_mod][
            (vec_log2fc_mod[this_mod] >= bin_start)
            *(vec_log2fc_mod[this_mod] < bin_end)
            ]
        bin_counts.append(len(binned_protein))
        binned_protein_log2fc.append(binned_protein)
        # binned_protein_log2fc.append(np.mean(binned_protein))
    mask = np.array([len(x) > 0 for x in binned_protein_log2fc])
    medians = [np.median(this_binned_log2fc) for this_binned_log2fc in binned_protein_log2fc]
    means = [np.mean(this_binned_log2fc) for this_binned_log2fc in binned_protein_log2fc]
    mod_means[this_mod] = means
    violin_parts = plt.violinplot(
        [this_bin for this_ind, this_bin in enumerate(binned_protein_log2fc) if mask[this_ind]],
        bin_centers[mask],
        # showmeans=True
    )
    for pc in violin_parts['bodies']:
        pc.set_facecolor(mod_color[this_mod])
        pc.set_edgecolor(mod_color[this_mod])
    violin_parts['cmaxes'].set_edgecolor(mod_color[this_mod])
    # violin_parts['cmeans'].set_edgecolor(mod_color[this_mod])
    violin_parts['cmins'].set_edgecolor(mod_color[this_mod])
    violin_parts['cbars'].set_edgecolor(mod_color[this_mod])
    plt.plot(bin_centers, means, c=mod_color[this_mod])
    plt.plot(bin_centers, means, c=mod_color[this_mod], marker='o')
    total_counts = np.sum(bin_counts)
    # plt.axhline(y=0, c='gray')
    # plt.bar(bin_centers, binned_protein_log2fc, width=bin_width*0.8)
    # for this_bin_center, this_bin_count, this_protein_log2fc in zip(bin_centers, bin_counts, binned_protein_log2fc):
    #     if not np.isnan(this_protein_log2fc):
    #         plt.text(this_bin_center, this_protein_log2fc + np.sign(this_protein_log2fc) * 0.025, this_bin_count,
    #              horizontalalignment='center', verticalalignment='center')
    plt.title(f'N={total_counts}')
    plt.xlim([-max_bin, max_bin])
    plt.xticks(bin_edges)
    # plt.ylim([-max_y, max_y])
    max_y = (np.nanmax(np.abs(np.concatenate(binned_protein_log2fc))) // 0.5 + 1) * 0.5
    plt.ylim([-max_y, max_y])
    plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
    # if mod_ind==0:
    plt.ylabel('$log_{2}FC (protein)$', fontsize=12)
plt.subplot(1, 3, 3)
for this_mod in mods:
    plt.plot(bin_centers, mod_means[this_mod], mod_color[this_mod], label=f'${dict_mod_display[this_mod]}$')
    plt.plot(bin_centers, mod_means[this_mod], mod_color[this_mod], marker='o')
plt.xlabel(f'$log_{2}FC (mod)$', fontsize=12)
plt.ylabel('$log_{2}FC (protein)$', fontsize=12)
ymax = 0.2
plt.ylim([-ymax, ymax])
plt.xticks(bin_edges)
plt.yticks(np.arange(-ymax, ymax+0.1, 0.1))
plt.legend(fontsize=12)
plt.suptitle(f'{this_day}', fontsize=15)
plt.savefig(os.path.join(img_out, f'log2fc_protein_vs_mods_{this_day}.png'), bbox_inches='tight')

# ### KS statistics ###
# gene_ks_stats = {}
# for gene, mod_ratios in gene_mod_values.items():
#         gene_ks_stats[gene] = {}
#         for this_mod in ['m6A', 'psi']:
#             if np.concatenate(list(mod_ratios[this_mod].values())).any():
#                 num_sites = len(mod_ratios[this_mod][f'SHAM_{this_day}'])
#                 sham_vals = mod_ratios[this_mod][f'SHAM_{this_day}']
#                 tac_vals = mod_ratios[this_mod][f'TAC_{this_day}']
#                 ks_stats = kstest(sham_vals, tac_vals)
#                 gene_ks_stats[gene][this_mod] = (num_sites, ks_stats.statistic, ks_stats.pvalue,
#                                                  np.mean(sham_vals), np.mean(tac_vals), np.mean(tac_vals)-np.mean(sham_vals))
#
#
# filtered_gene_ks_stats = {this_gene: {this_mod: this_ks for this_mod, this_ks in this_mod_vals.items()
#                                       if (this_ks[0]>=thresh_num_sites) and (this_ks[2]<thresh_pval)}
#                           for this_gene, this_mod_vals in gene_ks_stats.items()}
# filtered_gene_ks_stats = {k: v for k, v in filtered_gene_ks_stats.items() if len(v)}
# # for k, v in filtered_gene_ks_stats.items():
# #     print(k, v)
#
# df_ks_stats = pd.DataFrame(
#     [
#         [gene, this_mod, *this_mod_stats]
#         for gene, stats in filtered_gene_ks_stats.items() for this_mod, this_mod_stats in stats.items()
#     ],
#     columns=['gene', 'mod', 'num_sites', 'ks', 'pval', 'mean_SHAM', 'mean_TAC', 'delta']
# )
# df_ks_stats.sort_values('pval', inplace=True)
# df_ks_stats.to_csv(os.path.join(img_out, f'filtered_ks_stats_{this_day}.tsv'), sep='\t', float_format='%.3f', index=False)
#
# for this_gene, _ in filtered_gene_ks_stats.items():
#     this_ks_stats = gene_ks_stats[this_gene]
#     this_gene_mod_sites = gene_mod_values[this_gene]
#     plt.figure(figsize=(10, 4))
#     for mod_ind, this_mod in enumerate(mods):
#         plt.subplot(1, 2, mod_ind+1)
#         for cond_ind, this_cond in enumerate(conditions):
#             violin_parts = plt.violinplot(this_gene_mod_sites[this_mod][f'{this_cond}_{this_day}'],
#                                           [cond_ind], showmedians=True, widths=1)
#             for pc in violin_parts['bodies']:
#                 pc.set_facecolor(condition_color[this_cond])
#                 pc.set_edgecolor(condition_color[this_cond])
#             violin_parts['cmaxes'].set_edgecolor(condition_color[this_cond])
#             violin_parts['cmedians'].set_edgecolor(condition_color[this_cond])
#             violin_parts['cmins'].set_edgecolor(condition_color[this_cond])
#             violin_parts['cbars'].set_edgecolor(condition_color[this_cond])
#         plt.xticks(range(len(conditions)), conditions)
#         plt.ylim([0, 100])
#         plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
#         plt.title(f'ks={this_ks_stats[this_mod][1]:.2f}, pval={this_ks_stats[this_mod][2]:.2f}')
#     plt.suptitle(f'{this_gene}\n{this_day}', fontsize=15)
#     plt.savefig(os.path.join(img_out, f'{this_gene}_{this_day}.png'), bbox_inches='tight')
#     plt.close('all')
#
# ### <S> vs delta<S>, volcano plot ###
# gene_mod_delta = {}
# for gene, mod_values in gene_mod_values.items():
#     for this_mod, this_values in mod_values.items():
#         if len(this_values[f'SHAM_{this_day}'])>=thresh_num_sites:
#             mean_sham = np.mean(this_values[f'SHAM_{this_day}'])
#             mean_tac = np.mean(this_values[f'TAC_{this_day}'])
#             if gene not in gene_mod_delta.keys():
#                 gene_mod_delta[gene] = {}
#             gene_mod_delta[gene][this_mod] = (mean_sham, mean_tac-mean_sham)
#
# plt.figure(figsize=(10, 4))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     mean_s_sham, delta_mean_s = np.vstack(
#         [mod_delta[this_mod] for gene, mod_delta in gene_mod_delta.items() if mod_delta.get(this_mod) is not None]
#     ).T
#     plt.scatter(delta_mean_s, mean_s_sham, s=3)
#     plt.xlim([-20, 20])
#     plt.ylim([0, 100])
#     plt.ylabel(rf'$\langle S_{{{dict_mod_display[this_mod]}}} \rangle (SHAM)$', fontsize=12)
#     plt.xlabel(rf'$\langle S_{{{dict_mod_display[this_mod]}}} \rangle (TAC - SHAM)$', fontsize=12)
# plt.suptitle(f'{this_day}', fontsize=15)
# plt.savefig(os.path.join(img_out, f'volcano_{this_day}.png'), bbox_inches='tight')