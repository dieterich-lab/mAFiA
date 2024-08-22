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

# ds = 'TAC_SHAM_merged_TAC_merged',
# ds = 'Diet_WT_CD_merged_WT_WD_merged',
ds = 'HFpEF_ctrl_merged_HFpEF_merged'

ds_names = {
    'TAC_SHAM_merged_TAC_merged': 'TAC',
    'Diet_WT_CD_merged_WT_WD_merged': 'Diet',
    'HFpEF_ctrl_merged_HFpEF_merged': 'HFpEF',
}

ds_colors = {
    this_ds: this_color for (this_ds, this_color) in zip([
        'TAC_SHAM_merged_TAC_merged',
        'Diet_WT_CD_merged_WT_WD_merged',
        'HFpEF_ctrl_merged_HFpEF_merged',
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

df = pd.read_csv(os.path.join(ks_dir, f'{ds}.bed'), sep='\t', dtype={'chrom': str}).drop(columns=['score'])
df_thresh = df[df['pval'] < thresh_pval]
df_thresh['log2fc'] = np.log2(df_thresh['modRatio_2'] / df_thresh['modRatio_1'])

annotated_sites_file = os.path.join(img_out, f'annotated_diff_sites_{ds_names[ds]}.tsv')
if os.path.exists(annotated_sites_file):
    df_annotated_site = pd.read_csv(annotated_sites_file, sep='\t', dtype={'chrom': str})
else:
    df_annotated_sites = get_df_annotated(df_thresh)
    df_annotated_sites.to_csv(annotated_sites_file, sep='\t', index=False)

thresh_log2fc = 0.25
thresh_min_num_sites = 5
display_genes = 20
for pos_neg in ['positive', 'negative']:
    plt.figure(figsize=(15, 15))
    for mod_ind, this_mod in enumerate(mods):
        this_mod_df_annotated_sites = df_annotated_sites[df_annotated_sites['name']==this_mod]
        gene_ds_delta = {}
        for this_gene, site_counts in Counter(this_mod_df_annotated_sites['gene']).most_common():
            sub_df = this_mod_df_annotated_sites[this_mod_df_annotated_sites['gene'] == this_gene]
            if (len(sub_df) >= thresh_min_num_sites) and (np.abs(sub_df['log2fc'].values).max() >= thresh_log2fc):
                gene_ds_delta[this_gene] = [x for x in sub_df['log2fc'].values if ~np.isnan(x) and ~np.isinf(x)]

        gene_ds_delta = {k: v for k, v in gene_ds_delta.items() if len(v)}
        gene_ds_delta = {k: v
                         for k, v in sorted(gene_ds_delta.items(), key=lambda item: np.mean(item[1]), reverse=True)}

        if pos_neg == 'positive':
            sorted_genes = list(gene_ds_delta.keys())[:display_genes]
        elif pos_neg == 'negative':
            sorted_genes = list(gene_ds_delta.keys())[::-1][:display_genes]

        plt.subplot(2, 1, mod_ind+1)
        for col_ind, this_gene in enumerate(sorted_genes):
            ds_delta = gene_ds_delta[this_gene]
            plt.scatter([col_ind] * len(ds_delta), ds_delta, c=ds_colors[ds])
        plt.xticks(range(len(sorted_genes)), sorted_genes)
        plt.ylim([-4, 4])
        plt.axhline(y=0, c='gray', alpha=0.5, ls='--')
        plt.xlabel('Gene', fontsize=12)
        plt.ylabel(rf'$log2fc({{{dict_mod_display[this_mod]}}})$', fontsize=12)
        plt.title(f'${{{dict_mod_display[this_mod]}}}$', fontsize=15)
    plt.suptitle(f'{ds_names[ds]} \n {display_genes} genes with most {pos_neg} log2fc', fontsize=20)
    plt.savefig(os.path.join(img_out, f'individual_sites_{pos_neg}_log2fc_{ds_names[ds]}.png'), bbox_inches='tight')


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
