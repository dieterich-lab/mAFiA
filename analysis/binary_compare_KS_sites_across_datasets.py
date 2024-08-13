import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import numpy as np
from tqdm import tqdm
import pybedtools
from collections import Counter
from functools import reduce
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def get_df_annotated(in_df):
    df_annotated = []
    for _, this_row in tqdm(in_df.iterrows()):
        # this_intersect_df = annot.intersect(pybedtools.BedTool.from_dataframe(pd.DataFrame(this_row).T), s=True, wa=True).to_dataframe()
        this_intersect_df = annot.intersect(pybedtools.BedTool.from_dataframe(pd.DataFrame(this_row).T), wa=True).to_dataframe()
        if len(this_intersect_df)==0:
            continue
        this_intersect_df = this_intersect_df[this_intersect_df['strand'] == this_row['strand']]
        if len(this_intersect_df)==0:
            continue
        if len(this_intersect_df['name'].unique())>1:
            common_names =[k for k, v in Counter(this_intersect_df['name']).items() if v > 1]
            if len(common_names)>1:
                print('Non-unique gene:\n')
                print(this_intersect_df)
                continue
            else:
                this_intersect_df = this_intersect_df[this_intersect_df['name'] == common_names[0]]
        this_intersect_df.drop(columns=['itemRgb'], inplace=True)
        this_intersect_df.rename(columns={'thickStart': 'source', 'thickEnd': 'feature_type', 'blockCount': 'description'}, inplace=True)

        # print(this_row)
        # print(this_intersect_df)

        this_gene_id = this_intersect_df['name'].values[0]
        this_gene_name = [this_field.split(' ')[1]
                          for this_field in this_intersect_df[this_intersect_df['feature_type']=='gene']['description'].iloc[0].split('; ')
                          if this_field.split(' ')[0]=='gene_name'][0].strip('\"')
        this_gene_start, this_gene_end = this_intersect_df[this_intersect_df['feature_type']=='gene'][['start', 'end']].values[0]

        feature_type = [this_feature for this_feature in this_intersect_df['feature_type'].unique() if this_feature not in ['gene', 'transcript', 'exon']]
        if len(feature_type)==0:
            this_feature = None
        elif len(feature_type)==1:
            this_feature = feature_type[0]
        else:
            this_feature = ', '.join(feature_type)

        this_row['gene'] = this_gene_name
        this_row['geneStart'] = this_gene_start
        this_row['geneEnd'] = this_gene_end
        this_row['feature'] = this_feature

        df_annotated.append(this_row)

    df_annotated = pd.DataFrame(df_annotated)
    # df_annotated.to_csv(os.path.join(ks_dir, f'annotated_sites_{day}.tsv'), sep='\t', index=False)
    return df_annotated


mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}


ks_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/KS_test'
img_out = '/home/adrian/img_out/comparison_KS_sites_across_datasets'
os.makedirs(img_out, exist_ok=True)

annot_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.bed'
df_annot = pd.read_csv(annot_file, sep='\t', names=[
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'other',
    'description',
], dtype={'chrom': str})
df_annot = df_annot[df_annot['source']=='ensembl_havana']
annot = pybedtools.BedTool.from_dataframe(df_annot)

ds = [
    # 'TAC_SHAM_day7_TAC_day7',
    'Diet_WT_CD_merged_WT_WD_merged',
    'HFpEF_ctrl_merged_HFpEF_merged'
]

ds_names = {
    'TAC_SHAM_day7_TAC_day7': 'TAC',
    'Diet_WT_CD_merged_WT_WD_merged': 'Diet',
    'HFpEF_ctrl_merged_HFpEF_merged': 'HFpEF'
}

ds_colors = {
    this_ds: this_color for (this_ds, this_color) in zip([
        'TAC_SHAM_day7_TAC_day7',
        'Diet_WT_CD_merged_WT_WD_merged',
        'HFpEF_ctrl_merged_HFpEF_merged'
    ], ['g', 'b', 'r'])
}

merge_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'strand',
    'ref5mer'
]

rename_fields = [
    'ks_stat',
    'pval',
    'coverage_1',
    'coverage_2',
    'modRatio_1',
    'modRatio_2',
    'delta',
    'log2fc'
]

thresh_pval = 0.05

all_dfs = []
for this_ds in ds:
    this_df = pd.read_csv(os.path.join(ks_dir, f'{this_ds}.bed'), sep='\t', dtype={'chrom': str}).drop(columns=['score'])
    this_df_thresh = this_df[this_df['pval'] < thresh_pval]
    this_df_thresh['log2fc'] = np.log2(this_df_thresh['modRatio_2'] / this_df_thresh['modRatio_1'])
    this_df_thresh.rename(
        columns={this_field: this_field + f'_{this_ds}' for this_field in rename_fields},
        inplace=True
    )
    all_dfs.append(this_df_thresh)

df_merged = reduce(lambda left, right: pd.merge(left, right, on=merge_fields), all_dfs)
df_annotated_sites = get_df_annotated(df_merged)

thresh_log2fc = 0.25
thresh_min_num_sites = 2

plt.figure(figsize=(15, 15))
for mod_ind, this_mod in enumerate(mods):
    this_mod_df_annotated_sites = df_annotated_sites[df_annotated_sites['name']==this_mod]
    gene_ds_delta = {}
    for this_gene, site_counts in Counter(this_mod_df_annotated_sites['gene']).most_common():
        sub_df = this_mod_df_annotated_sites[this_mod_df_annotated_sites['gene'] == this_gene]
        if (len(sub_df) >= thresh_min_num_sites) \
                and (np.abs(sub_df[[f'log2fc_{ds[0]}', f'log2fc_{ds[1]}']].values).max() >= thresh_log2fc):
            gene_ds_delta[this_gene] = [sub_df[f'log2fc_{ds[0]}'].values, sub_df[f'log2fc_{ds[1]}'].values]

    # sorted_genes = np.array(list(gene_ds_delta.keys()))[
    #     np.argsort([(np.mean(this_val[1]) - np.mean(this_val[0])) for this_val in gene_ds_delta.values()])
    # ][::-1]
    sorted_genes = np.array(list(gene_ds_delta.keys()))[
        np.argsort([np.mean(this_val[0]) for this_val in gene_ds_delta.values()])
    ]

    plt.subplot(1, 2, mod_ind+1)
    labelled = False
    for row_ind, this_gene in enumerate(sorted_genes):
        ds_delta = gene_ds_delta[this_gene]
        if not labelled:
            for ds_ind, (this_ds, this_ds_delta) in enumerate(zip(ds, ds_delta)):
                plt.scatter(this_ds_delta, [row_ind-0.1+0.2*ds_ind] * len(this_ds_delta),
                            c=ds_colors[this_ds], label=ds_names[this_ds])
                labelled = True
        else:
            for ds_ind, (this_ds, this_ds_delta) in enumerate(zip(ds, ds_delta)):
                plt.scatter(this_ds_delta, [row_ind-0.1+0.2*ds_ind] * len(this_ds_delta),
                            c=ds_colors[this_ds])
    plt.gca().invert_yaxis()
    plt.yticks(range(len(gene_ds_delta.keys())), gene_ds_delta.keys())
    plt.xlim([-1.5, 1.5])
    plt.axvline(x=0, c='gray', alpha=0.5, ls='--')
    plt.legend(loc='lower left')
    plt.xlabel(rf'$log2fc$', fontsize=12)
    plt.ylabel('Gene', fontsize=12)
    plt.title(f'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
plt.suptitle(rf'$log2fc\geq{thresh_log2fc}$' + f'\nmin. {thresh_min_num_sites} sites per gene', fontsize=20)
plt.savefig(os.path.join(img_out, f'individual_sites_{ds_names[ds[0]]}_{ds_names[ds[1]]}.png'), bbox_inches='tight')

### meta-transcript profile ###
ds_hist = {
    this_ds: {
        this_mod: [] for this_mod in mods
    } for this_ds in ds
}
for this_gene in df_annotated_sites['gene'].unique():
    this_gene_df = df_annotated_sites[df_annotated_sites['gene'] == this_gene]
    this_gene_start = this_gene_df['geneStart'].values[0]
    this_gene_end = this_gene_df['geneEnd'].values[0]
    for this_mod in mods:
        this_gene_mod_df = this_gene_df[this_gene_df['name'] == this_mod]
        if len(this_gene_mod_df)==0:
            continue
        vec_pos = this_gene_mod_df['chromStart'].values
        if this_gene_mod_df['strand'].unique()[0]=='+':
            vec_norm_pos = (vec_pos-this_gene_start) / (this_gene_end-this_gene_start)
        else:
            vec_norm_pos = (this_gene_end-vec_pos) / (this_gene_end-this_gene_start)
        for this_ds in ds:
            ds_hist[this_ds][this_mod].extend(list(zip(vec_norm_pos, this_gene_df[f'log2fc_{this_ds}'])))

num_bins = 10
bin_edges = np.linspace(0, 1, num_bins+1)
bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5

ylim = [-1.5, 1.5]
plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    for this_ds in ds:
        vec_pos, vec_log2fc = np.vstack(ds_hist[this_ds][this_mod]).T
        vec_hist = []
        for bin_ind in range(len(bin_edges)-1):
            bin_start = bin_edges[bin_ind]
            bin_end = bin_edges[bin_ind+1]
            mask = (vec_pos >= bin_start) * (vec_pos<bin_end)
            masked_vals = [x for x in vec_log2fc[mask] if ~np.isnan(x) and ~np.isinf(x)]
            vec_hist.append(np.mean(masked_vals))
        plt.bar(bin_centers, vec_hist, width=0.08, label=ds_names[this_ds], color=ds_colors[this_ds], alpha=0.5)
    plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
    plt.ylim(ylim)
    plt.xlabel('% gene span', fontsize=12)
    plt.ylabel('Avg. log2fc', fontsize=12)
    plt.legend(loc='lower right')
plt.savefig(os.path.join(img_out, f'metagene_profile_{ds_names[ds[0]]}_{ds_names[ds[1]]}.png'), bbox_inches='tight')

### violinplot by transcript region ###
regions = [
    'five_prime_utr',
    'CDS',
    'three_prime_utr'
]

ylim = [-2, 2]

plt.figure(figsize=(10, 4))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    plt.axhline(y=0, c='gray', ls='--', alpha=0.2)
    labels = [(mpatches.Patch(color=ds_colors[this_ds]), ds_names[this_ds]) for this_ds in ds]
    for ds_ind, this_ds in enumerate(ds):
        region_log2fc = [
            df_annotated_sites[
                (df_annotated_sites['name'] == this_mod)
                * (df_annotated_sites['feature'] == this_region)
            ][f'log2fc_{this_ds}'].values
            for this_region in regions
        ]
        region_log2fc = [
            [x for x in this_region_log2fc if ~np.isnan(x) and ~np.isinf(x)]
            for this_region_log2fc in region_log2fc
        ]
        # region_num_sites = [len(x) for x in region_log2fc]
        # region_max = [max(x) for x in region_log2fc]
        # region_log2fc = [x if len(x) else [0] for x in region_log2fc]

        bplot = plt.boxplot(region_log2fc, positions=np.arange(len(regions))-0.2+ds_ind*0.4,
                    showfliers=False, patch_artist=True)
        for patch in bplot['boxes']:
            patch.set_facecolor(ds_colors[this_ds])
        for patch in bplot['medians']:
            patch.set_color('k')

        # violin_parts = plt.violinplot(region_log2fc, np.arange(len(regions))-0.2+ds_ind*0.4,
        #                               widths=0.3, showmedians=True)
        # for pc in violin_parts['bodies']:
        #     pc.set_facecolor(ds_colors[this_ds])
        #     pc.set_edgecolor(ds_colors[this_ds])
        # violin_parts['cmaxes'].set_edgecolor(ds_colors[this_ds])
        # violin_parts['cmedians'].set_edgecolor(ds_colors[this_ds])
        # violin_parts['cmins'].set_edgecolor(ds_colors[this_ds])
        # violin_parts['cbars'].set_edgecolor(ds_colors[this_ds])

        # for region_ind, (this_max, this_num_sites) in enumerate(zip(region_max, region_num_sites)):
        #     plt.text(region_ind-0.2+ds_ind*0.4, this_max+0.1, this_num_sites,
        #              c=ds_colors[this_ds], horizontalalignment='center')

    plt.legend(*zip(*labels), loc='upper left')
    plt.xticks(range(len(regions)), regions)
    plt.ylim(ylim)
    plt.xlabel('Region', fontsize=12)
    plt.ylabel('log2fc', fontsize=12)
    plt.title(rf'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
plt.savefig(os.path.join(img_out, f'boxplot_by_region_{ds_names[ds[0]]}_{ds_names[ds[1]]}.png'), bbox_inches='tight')
