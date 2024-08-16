import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import pybedtools
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm

day = 'day7'
thresh_pval = 0.05

ks_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/polyA/'
ks_file = os.path.join(ks_dir, f'KSTest_SHAM_{day}_below90_above90.bed')
img_out = '/home/adrian/img_out/single_site_ks_test'
df_ks = pd.read_csv(ks_file, sep='\t')
df_ks_sel = df_ks[(df_ks['pval']<thresh_pval)]
# ks_sel = pybedtools.BedTool.from_dataframe(df_ks_sel)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

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

def get_df_annotated():
    df_annotated = []
    for _, this_row in tqdm(df_ks_sel.iterrows()):
        this_intersect_df = annot.intersect(pybedtools.BedTool.from_dataframe(pd.DataFrame(this_row).T), s=True, wa=True).to_dataframe()
        if len(this_intersect_df)==0:
            continue
        if len(this_intersect_df['name'].unique())>1:
            continue
        this_intersect_df.drop(columns=['itemRgb'], inplace=True)
        this_intersect_df.rename(columns={'thickStart': 'source', 'thickEnd': 'feature_type', 'blockCount': 'description'}, inplace=True)

        # print(this_row)
        # print(this_intersect_df)

        this_gene_id = this_intersect_df['name'][0]
        this_gene_name = [this_field.split(' ')[1]
                          for this_field in this_intersect_df[this_intersect_df['feature_type']=='gene']['description'].iloc[0].split('; ')
                          if this_field.split(' ')[0]=='gene_name'][0].strip('\"')
        this_gene_start, this_gene_end = this_intersect_df[this_intersect_df['feature_type']=='gene'][['start', 'end']].values[0]

        feature_type = [this_feature for this_feature in this_intersect_df['feature_type'].unique() if this_feature not in ['gene', 'transcript', 'exon']]
        if len(feature_type)==0:
            this_feature = None
        elif len(feature_type)==1:
            this_feature = feature_type[0]
            # feature_start, feature_end = this_intersect_df[this_intersect_df['feature_type']==this_feature][['start', 'end']].values[0]
        else:
            this_feature = ', '.join(feature_type)

        # tx_ids = [this_field.split(' ')[1].strip('\"')
        #     for this_desc in this_intersect_df[this_intersect_df['feature_type'] == 'transcript']['description']
        #     for this_field in this_desc.split('; ')
        #     if this_field.split(' ')[0] == 'transcript_id'
        # ]
        # print(tx_ids)

        this_row['gene'] = this_gene_name
        this_row['geneStart'] = this_gene_start
        this_row['geneEnd'] = this_gene_end
        this_row['feature'] = this_feature

        df_annotated.append(this_row)

        # tx_regions = []
        # for _, this_annot in this_intersect_df.iterrows():
        #     if this_annot['feature_type'] not in ['gene', 'transcript', 'exon']:
        #         this_out_row = this_row.copy()
        #
        #         this_tx_id = [this_field.split(' ')[1].strip('\"')
        #                       for this_field in this_annot['description'].split('; ')
        #                       if this_field.split(' ')[0] == 'transcript_id'][0]
        #         # if len(this_tx_id)>1:
        #         #     print('>1 tx_ids!')
        #         # this_tx_id = this_tx_id[0]
        #         this_feature_start, this_feature_end = this_annot[['start', 'end']].values
        #         this_out_row['transcript'] = this_tx_id
        #         this_out_row['feature'] = this_annot['feature_type']
        #         this_out_row['featureStart'] = this_feature_start
        #         this_out_row['featureEnd'] = this_feature_end
        #
        #         df_annotated.append(this_out_row)

    df_annotated = pd.DataFrame(df_annotated)
    df_annotated.to_csv(os.path.join(ks_dir, f'annotated_sites_{day}.tsv'), sep='\t', index=False)
    return df_annotated

annotated_sites_file = os.path.join(ks_dir, f'annotated_sites_{day}.tsv')
if os.path.exists(annotated_sites_file):
    df_annotated = pd.read_csv(annotated_sites_file, sep='\t')
else:
    df_annotated = get_df_annotated()


regions = ['five_prime_utr', 'CDS', 'stop_codon', 'three_prime_utr']

plt.figure(figsize=(12, 5))
for this_mod_ind, this_mod in enumerate(mods):
    df_mod = df_annotated[df_annotated['name'] == this_mod]
    region_mod_vals = [df_mod[df_mod['feature'] == this_region]['delta'].values for this_region in regions]
    region_mod_vals = [this_region_mod_vals if len(this_region_mod_vals) else [0] for this_region_mod_vals in region_mod_vals]

    plt.subplot(1, 2, this_mod_ind+1)
    plt.violinplot(region_mod_vals, range(len(regions)), showmeans=True)
    plt.ylim([-25, 25])
    plt.xlabel('Region', fontsize=12)
    plt.ylabel(f'$\Delta S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
    plt.axhline(y=0, color='black')
    plt.xticks(range(len(regions)), regions)
plt.suptitle(f'TAC vs. SHAM\n{day}', fontsize=15)
plt.savefig(os.path.join(img_out, 'mod_dist_by_region.png'), bbox_inches='tight')


### calculate gene position ###
# num_bins = {
#     'five_prime_utr': 2,
#     'CDS': 5,
#     'stop_codon': 1,
#     'three_prime_utr': 2
# }
# total_num_bins = np.sum(list(num_bins.values()))
# bin_centers = np.arange(total_num_bins) + 0.5
# bin_width = bin_centers[1] - bin_centers[0]
# features = num_bins.keys()
#
# mod_feature_pos_delta = {this_mod: {this_feature: [] for this_feature in features} for this_mod in mods}
# # for _, this_row in df_annotated.iterrows():
# for this_chromStart in df_annotated['chromStart'].unique():
#     sub_df = df_annotated[df_annotated['chromStart']==this_chromStart]
#     if len(sub_df)>1:
#         if (len(sub_df['feature'].unique()) == 1) \
#                 and (len(sub_df['featureStart'].unique()) == 1) \
#                 and (len(sub_df['featureEnd'].unique()) == 1):
#             # print(sub_df)
#             this_feature_start = sub_df['featureStart'].values[0]
#             this_feature_end = sub_df['featureEnd'].values[0]
#             feature_len = this_feature_end - this_feature_start
#             this_strand = sub_df['strand'].unique()[0]
#             if this_strand =='+':
#                 site_dist = this_chromStart - this_feature_start
#             else:
#                 site_dist = this_feature_end - this_chromStart + 1
#             site_dist_norm = site_dist / feature_len
#             mod_feature_pos_delta[sub_df['name'].values[0]][sub_df['feature'].unique()[0]].append(
#                 (site_dist_norm, sub_df['delta'].values[0])
#             )
#
# feature_bin_edges = {}
# for this_feature in features:
#     feature_bin_edges[this_feature] = np.linspace(0, 1, num_bins[this_feature]+1)
#
# ylim = [-10, 10]
#
# plt.figure(figsize=(10, 4))
# for this_mod_ind, this_mod in enumerate(mods):
#     all_feature_binned_deltas = []
#     for this_feature_ind, this_feature in enumerate(features):
#         if len(mod_feature_pos_delta[this_mod][this_feature])==0:
#             binned_deltas = [0] * (len(feature_bin_edges[this_feature]) - 1)
#         else:
#             vec_pos, vec_delta = np.vstack(mod_feature_pos_delta[this_mod][this_feature]).T
#             binned_deltas = []
#             bin_edges = feature_bin_edges[this_feature]
#             for bin_i in range(len(bin_edges)-1):
#                 bin_start = bin_edges[bin_i]
#                 bin_end = bin_edges[bin_i+1]
#                 mask = (vec_pos>=bin_start) * (vec_pos<bin_end)
#                 binned_deltas.append(np.mean(vec_delta[mask]))
#         all_feature_binned_deltas.extend(binned_deltas)
#     plt.subplot(1, 2, this_mod_ind+1)
#     plt.bar(bin_centers, all_feature_binned_deltas, width=bin_width*0.8)
#     plt.ylim(ylim)
#     plt.xlabel('Metagene position', fontsize=12)
#     plt.ylabel(f'$\Delta S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
#     plt.axvspan(
#         num_bins['five_prime_utr'],
#         num_bins['five_prime_utr'] + num_bins['CDS'] + num_bins['stop_codon'],
#         color='gray', alpha=0.2)
#     plt.axhline(y=0, color='black')
#     plt.xticks([num_bins['five_prime_utr']/2, num_bins['five_prime_utr']+(num_bins['CDS'] + num_bins['stop_codon'])/2, num_bins['five_prime_utr'] + num_bins['CDS'] + num_bins['stop_codon'] + num_bins['three_prime_utr']/2], ["5' UTR", "CDS", "3' UTR"])
# plt.suptitle(day, fontsize=15)
# plt.savefig(os.path.join(img_out, f'metagene_profile_{day}.png'), bbox_inches='tight')
