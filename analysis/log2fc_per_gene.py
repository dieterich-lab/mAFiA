import pandas as pd
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import kstest
import numpy as np
from scipy.stats import trim_mean

this_day = 'day21'

thresh_confidence = 0.0
thresh_coverage = 20
thresh_pval = 0.1
thresh_num_sites = 10
mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
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

proteome_xls = '/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis.xlsx'
df_proteome = pd.read_excel('/home/adrian/Data/Dewenter_TAC_Backs_lab/TACOMA_differential_proteome_analysis.xlsx')
gene_protein_log2fc = {}
for this_gene in df_proteome['SYMBOL'].unique():
    sub_df = df_proteome[
        (df_proteome['SYMBOL']==this_gene)
        * (df_proteome['main.comparison']=='TAC vs Sham')
        * (df_proteome['tp']==('Day ' + this_day.lstrip('day')))
    ]
    gene_protein_log2fc[this_gene] = sub_df['logFC'].values[0]

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'
conditions = ['SHAM', 'TAC']
condition_days = [f'{this_cond}_{this_day}' for this_cond in conditions]

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

### segregate by gene ###
gene_mod_values = {}
for this_gene in tqdm(df_gene_bed['name'].unique()):
    gene_subdf = df_gene_bed[df_gene_bed['name']==this_gene]
    chrom = gene_subdf['chrom'].unique()[0]
    chromStart = gene_subdf['chromStart'].min()
    chromEnd = gene_subdf['chromEnd'].max()
    strand = gene_subdf['strand'].unique()[0]
    df_gene_mod_sites = df_merged[
        (df_merged['chrom']==chrom)
        * (df_merged['chromStart']>=chromStart)
        * (df_merged['chromEnd']<chromEnd)
        * (df_merged['strand']==strand)
    ]
    this_gene_mod_values = {
        this_mod: {
            this_cond_day: df_gene_mod_sites[df_gene_mod_sites['name']==this_mod][f'modRatio_{this_cond_day}'].values
            for this_cond_day in condition_days
        }
        for this_mod in ['m6A', 'psi']
    }

    if np.concatenate([this_cond_vals for this_mod_vals in this_gene_mod_values.values() for this_cond_vals in this_mod_vals.values()]).any():
        gene_mod_values[this_gene] = this_gene_mod_values


### log2FC ###
metric = 'mean'
gene_log2fc_stats = {}
for gene, mod_ratios in gene_mod_values.items():
        gene_log2fc_stats[gene] = {}
        for this_mod in ['m6A', 'psi']:
            if np.concatenate(list(mod_ratios[this_mod].values())).any():
                num_sites = len(mod_ratios[this_mod][f'SHAM_{this_day}'])
                sham_vals = mod_ratios[this_mod][f'SHAM_{this_day}']
                tac_vals = mod_ratios[this_mod][f'TAC_{this_day}']
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
                gene_log2fc_stats[gene][this_mod] = (num_sites, sham_metric, tac_metric, log2fc)

### all ###
max_range = 1.5
gene_log2fc_m6A = {this_gene: this_gene_log2fc.get('m6A')[3] for this_gene, this_gene_log2fc in gene_log2fc_stats.items() if this_gene_log2fc.get('m6A') is not None}
gene_log2fc_psi = {this_gene: this_gene_log2fc.get('psi')[3] for this_gene, this_gene_log2fc in gene_log2fc_stats.items() if this_gene_log2fc.get('psi') is not None}
all_log2fc = {
    'm6A': gene_log2fc_m6A.values(),
    'psi': gene_log2fc_psi.values()
}
plt.figure(figsize=(10, 5))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    plt.hist(all_log2fc[this_mod], bins=50, range=[-max_range, max_range], alpha=0.5)
    plt.axvline(x=0, c='g')
    plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
    plt.ylabel('Gene counts', fontsize=12)
    plt.title(f'N={len(all_log2fc[this_mod])}', fontsize=15)

### common gene and thresholded by min. sites ###
common_genes, common_log2fc_m6A, common_log2fc_psi = np.vstack(
    [(gene, val['m6A'][3], val['psi'][3]) for gene, val in gene_log2fc_stats.items()
     if (set(val.keys()).issuperset((set(['m6A', 'psi']))))
     and (val['m6A'][0]>=thresh_num_sites) and (val['psi'][0]>=thresh_num_sites)
     ]
).T
common_log2fc = {
    'm6A': np.float64(common_log2fc_m6A),
    'psi': np.float64(common_log2fc_psi)
}

plt.figure(figsize=(10, 5))
for mod_ind, this_mod in enumerate(mods):
    that_mod = [m for m in mods if m!=this_mod][0]
    plt.subplot(1, 2, mod_ind+1)
    plt.hist(common_log2fc[this_mod][common_log2fc[that_mod]>=0], bins=50, range=[-max_range, max_range], fc='b', label=f'${dict_mod_display[that_mod]}$ up', alpha=0.5)
    plt.hist(common_log2fc[this_mod][common_log2fc[that_mod]<0], bins=50, range=[-max_range, max_range], fc='r', label=f'${dict_mod_display[that_mod]}$ down', alpha=0.5)
    plt.axvline(x=0, c='g')
    plt.legend(loc='upper left')
    plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
    plt.ylabel('Gene counts', fontsize=12)
    plt.title(f'N={len(common_log2fc[this_mod])}', fontsize=15)


### protein vs mods ###
common_genes_m6A = list(set(gene_log2fc_m6A.keys()).intersection(set(gene_protein_log2fc.keys())))
vec_log2fc_m6A, vec_log2fc_m6A_protein = np.vstack(
    [(gene_log2fc_m6A[this_gene], gene_protein_log2fc[this_gene]) for this_gene in common_genes_m6A]
).T
common_genes_psi = list(set(gene_log2fc_psi.keys()).intersection(set(gene_protein_log2fc.keys())))
vec_log2fc_psi, vec_log2fc_psi_protein = np.vstack(
    [(gene_log2fc_psi[this_gene], gene_protein_log2fc[this_gene]) for this_gene in common_genes_psi]
).T
vec_log2fc_mod = {}
vec_log2fc_mod['m6A'] = vec_log2fc_m6A
vec_log2fc_mod['psi'] = vec_log2fc_psi
vec_log2fc_protein = {}
vec_log2fc_protein['m6A'] = vec_log2fc_m6A_protein
vec_log2fc_protein['psi'] = vec_log2fc_psi_protein

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    plt.scatter(vec_log2fc_mod[this_mod], vec_log2fc_protein[this_mod], s=2)
    plt.xlim([-3, 3])
    plt.ylim([-2, 2])
    plt.axvline(x=0, c='r', ls='--')
    plt.axhline(y=0, c='r', ls='--')
    plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
    plt.ylabel('$log_{2}FC (protein)$', fontsize=12)
plt.suptitle(f'{this_day}', fontsize=15)

### aggregate protein changes binned by mod changes ###
bins = np.linspace(-3, 3, 7)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    binned_protein_log2fc = []
    total_n = 0
    for i in range(len(bins)-1):
        bin_start = bins[i]
        bin_end = bins[i+1]
        binned_protein = vec_log2fc_protein[this_mod][
            (vec_log2fc_mod[this_mod] >= bin_start)
            *(vec_log2fc_mod[this_mod] < bin_end)
            ]
        total_n += len(binned_protein)
        binned_protein_log2fc.append(np.mean(binned_protein))
    plt.axhline(y=0, c='gray')
    plt.bar(bin_centers, binned_protein_log2fc)
    plt.title(f'N={total_n}')
    plt.xlim([-3, 3])
    plt.ylim([-0.25, 0.25])
    plt.xlabel(f'$log_{2}FC ({dict_mod_display[this_mod]})$', fontsize=12)
    if mod_ind==0:
        plt.ylabel('$log_{2}FC (protein)$', fontsize=12)
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