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
from matplotlib_venn import venn3


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
    'TAC_SHAM_day7_TAC_day7',
    'Diet_WT_CD_merged_WT_WD_merged',
    'HFpEF_ctrl_merged_HFpEF_merged'
]

ds_names = {
    'TAC_SHAM_day7_TAC_day7': 'TAC',
    'Diet_WT_CD_merged_WT_WD_merged': 'Diet',
    'HFpEF_ctrl_merged_HFpEF_merged': 'HFpEF'
}

ds_colors = {
    this_ds: this_color for (this_ds, this_color) in zip(ds, ['g', 'b', 'r'])
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

ds_sites = [
    set([tuple(this_val) for this_val in this_df[['chrom', 'chromStart']].values])
    for this_df in all_dfs
]

plt.figure(figsize=(5, 5))
venn3(ds_sites, ds_names.values(), set_colors=ds_colors.values())
plt.title('Common differential sites across datasets', fontsize=15)
plt.savefig(os.path.join(img_out, 'venn_diagram_all_ds.png'), bbox_inches='tight')

df_merged = reduce(lambda left, right: pd.merge(left, right, on=merge_fields, how='inner'), all_dfs)
df_annotated_sites = get_df_annotated(df_merged)

gene_counts = Counter(df_annotated_sites['gene']).most_common()
with open(os.path.join(img_out, 'common_site_genes_all_ds.tsv'), 'w') as f_out:
    for this_count in gene_counts:
        f_out.write(f'{this_count[0]}\t{this_count[1]}\n')