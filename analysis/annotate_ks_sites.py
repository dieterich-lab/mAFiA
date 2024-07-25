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

day = 'day1'
thresh_pval = 0.05

ks_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/KS_test'
ks_file = os.path.join(ks_dir, f'SHAM_TAC_{day}.bed')
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

        this_gene_id = this_intersect_df['name'][0]
        feature_type = [this_feature for this_feature in this_intersect_df['feature_type'].unique() if this_feature not in ['gene', 'transcript', 'exon']]
        if len(feature_type)==1:
            this_feature = feature_type[0]
            feature_start, feature_end = this_intersect_df[this_intersect_df['feature_type']==this_feature][['start', 'end']].values[0]
        else:
            # print(feature_type)
            continue

        this_gene_name = [this_field.split(' ')[1]
                          for this_field in this_intersect_df[this_intersect_df['feature_type']=='gene']['description'].iloc[0].split('; ')
                          if this_field.split(' ')[0]=='gene_name'][0].strip('\"')

        this_row['gene_name'] = this_gene_name
        this_row['feature'] = this_feature
        this_row['featureStart'] = feature_start
        this_row['featureEnd'] = feature_end
        df_annotated.append(this_row)
    df_annotated = pd.DataFrame(df_annotated)
    df_annotated.to_csv(os.path.join(ks_dir, f'annotated_sites_{day}.tsv'), sep='\t', index=False)
    return df_annotated

annotated_sites_file = os.path.join(ks_dir, f'annotated_sites_{day}.tsv')
if os.path.exists(annotated_sites_file):
    df_annotated = pd.read_csv(annotated_sites_file, sep='\t')
else:
    df_annotated = get_df_annotated()

### calculate gene position ###
num_bins = {
    'five_prime_utr': 2,
    'CDS': 5,
    'stop_codon': 1,
    'three_prime_utr': 2
}
total_num_bins = np.sum(list(num_bins.values()))
bin_centers = np.arange(total_num_bins) + 0.5
bin_width = bin_centers[1] - bin_centers[0]
features = num_bins.keys()

mod_feature_pos_delta = {this_mod: {this_feature: [] for this_feature in features} for this_mod in mods}
for _, this_row in df_annotated.iterrows():
    feature_len = this_row['featureEnd'] - this_row['featureStart']
    if this_row['strand']=='+':
        site_dist = this_row['chromStart'] - this_row['featureStart']
    else:
        site_dist = this_row['featureEnd'] - this_row['chromStart'] + 1
    site_dist_norm = site_dist / feature_len
    mod_feature_pos_delta[this_row['name']][this_row['feature']].append((site_dist_norm, this_row['delta']))

feature_bin_edges = {}
for this_feature in features:
    feature_bin_edges[this_feature] = np.linspace(0, 1, num_bins[this_feature]+1)

ylim = [-8, 8]

plt.figure(figsize=(10, 4))
for this_mod_ind, this_mod in enumerate(mods):
    all_feature_binned_deltas = []
    for this_feature_ind, this_feature in enumerate(features):
        vec_pos, vec_delta = np.vstack(mod_feature_pos_delta[this_mod][this_feature]).T
        binned_deltas = []
        bin_edges = feature_bin_edges[this_feature]
        for bin_i in range(len(bin_edges)-1):
            bin_start = bin_edges[bin_i]
            bin_end = bin_edges[bin_i+1]
            mask = (vec_pos>=bin_start) * (vec_pos<bin_end)
            binned_deltas.append(np.mean(vec_delta[mask]))
        all_feature_binned_deltas.extend(binned_deltas)
    plt.subplot(1, 2, this_mod_ind+1)
    plt.bar(bin_centers, all_feature_binned_deltas, width=bin_width*0.8)
    plt.ylim(ylim)
    plt.xlabel('Metagene position', fontsize=12)
    plt.ylabel(f'$\Delta S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
    plt.axvspan(
        num_bins['five_prime_utr'],
        num_bins['five_prime_utr'] + num_bins['CDS'] + num_bins['stop_codon'],
        color='gray', alpha=0.2)
    plt.axhline(y=0, color='black')
    plt.xticks([num_bins['five_prime_utr']/2, num_bins['five_prime_utr']+(num_bins['CDS'] + num_bins['stop_codon'])/2, num_bins['five_prime_utr'] + num_bins['CDS'] + num_bins['stop_codon'] + num_bins['three_prime_utr']/2], ["5' UTR", "CDS", "3' UTR"])
plt.suptitle(day, fontsize=15)
plt.savefig(os.path.join(img_out, f'metagene_profile_{day}.png'), bbox_inches='tight')
