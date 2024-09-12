import os
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm

irfinder_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/IRFinder'
mod_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

img_out = '/home/adrian/img_out/IR_junctions/Diet'
os.makedirs(img_out, exist_ok=True)

mods = ['m6A', 'psi']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

# ds = 'HFpEF'
# conditions = ['ctrl_merged', 'HFpEF_merged']
ds = 'Diet'
conditions = ['WT_CD_merged', 'WT_WD_merged']
# ds = 'TAC'
# conditions = ['SHAM_merged', 'TAC_merged']
cond_names = [this_cond.rstrip('_merged') for this_cond in conditions]

### delta IR ratio ###
# this_cond = conditions[0]
cond_df_irf = {}
for this_cond in conditions:
    df_irf = pd.read_csv(os.path.join(irfinder_dir, f'{ds}_{this_cond}.txt'), sep='\t', dtype={'Chr': str})
    df_irf = df_irf.rename(columns={
        'Chr': 'chrom',
        'Start': 'chromStart',
        'End': 'chromEnd',
        'Name': 'name',
        'Null': 'score',
        'Strand': 'strand'
    })
    df_irf_clean = df_irf.loc[['known-exon' not in this_name.split('/')[-1] for this_name in df_irf['name']], :]
    cond_df_irf[this_cond] = df_irf_clean[df_irf_clean['Warnings'] == '-']

thresh_IRratio = 0.05
thresh_IntronDepth = 5

df_irf_merged = pd.merge(*list(cond_df_irf.values()), on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'],
                         suffixes=['_'+this_cond_name for this_cond_name in cond_names])
df_irf_merged_thresh = df_irf_merged[
    ((df_irf_merged.loc[:, df_irf_merged.keys().str.find('IRratio') >= 0] >= thresh_IRratio).any(axis=1))
    * ((df_irf_merged.loc[:, df_irf_merged.keys().str.find('IntronDepth_') >= 0] >= thresh_IntronDepth).any(axis=1))
].copy()
df_irf_merged_thresh['delta_IRratio'] = df_irf_merged_thresh[f'IRratio_{cond_names[1]}'] - \
                                        df_irf_merged_thresh[f'IRratio_{cond_names[0]}']
df_irf_merged_thresh.sort_values('delta_IRratio', inplace=True)

df_irf_merged_thresh.to_csv(
    os.path.join(img_out, f"ir_list_{'_'.join(cond_names)}_IRratio{thresh_IRratio}_IntronDepth{thresh_IntronDepth}.tsv"),
    sep='\t',
    index=False
)

### scatter plot of IR ratios ###
xylim = [-0.05, 1.05]
plt.figure(figsize=(5, 5))
plt.plot(df_irf_merged_thresh[f'IRratio_{cond_names[0]}'], df_irf_merged_thresh[f'IRratio_{cond_names[1]}'], '.')
plt.plot([0, 1], [0, 1], c='gray', ls='--')
plt.xlim(xylim)
plt.ylim(xylim)
plt.xlabel(f'IR ratio, {cond_names[0]}')
plt.ylabel(f'IR ratio, {cond_names[1]}')

### delta mod ###
thresh_conf = 0.0
cond_df_mod = {}
for this_cond in conditions:
    cond_df_mod[this_cond] = pd.read_csv(os.path.join(mod_dir, ds, this_cond, f'chrALL.mAFiA.sites.bed'),
                                         sep='\t', dtype={'chrom': str})
df_mod_merged = pd.merge(*list(cond_df_mod.values()),  on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'],
                         suffixes=['_'+this_cond_name for this_cond_name in cond_names])
df_mod_merged_thresh = df_mod_merged[
    (df_mod_merged[f'confidence_{cond_names[0]}'] >= thresh_conf)
    * (df_mod_merged[f'confidence_{cond_names[1]}'] >= thresh_conf)
].copy()
df_mod_merged_thresh['delta_modRatio'] = df_mod_merged_thresh[f'modRatio_{cond_names[1]}'] - df_mod_merged_thresh[f'modRatio_{cond_names[0]}']


margin = 50

delta_ir_avg_mod_deltas = []
for _, this_row in tqdm(df_irf_merged_thresh.iterrows()):
    regions = [
        (this_row['chromStart']-margin, this_row['chromStart']),
        (this_row['chromEnd'], this_row['chromEnd']+margin)
        ]
    mod_deltas = {this_mod: [] for this_mod in mods}
    for this_region in regions:
        sub_df_mod = df_mod_merged_thresh[
            (df_mod_merged_thresh['chrom'] == this_row['chrom'])
            * (df_mod_merged_thresh['chromStart'] >= this_region[0])
            * (df_mod_merged_thresh['chromStart'] < this_region[1])
            * (df_mod_merged_thresh['strand'] == this_row['strand'])
        ]
        if len(sub_df_mod):
            for this_mod in mods:
                this_sub_df_mod = sub_df_mod[sub_df_mod['name'] == this_mod]
                if len(this_sub_df_mod):
                    mod_deltas[this_mod].extend(this_sub_df_mod['delta_modRatio'].values)
    # if (len(mod_deltas['m6A'])>=5):
    #     mean_mod_delta_m6A = np.mean(mod_deltas['m6A'])
    # else:
    #     mean_mod_delta_m6A = np.nan
    # if (len(mod_deltas['psi'])>=5):
    #     mean_mod_delta_psi = np.mean(mod_deltas['psi'])
    # else:
    #     mean_mod_delta_psi = np.nan
    mean_mod_delta_m6A = np.mean(mod_deltas['m6A'])
    mean_mod_delta_psi = np.mean(mod_deltas['psi'])
    delta_ir_avg_mod_deltas.append((this_row['delta_IRratio'], mean_mod_delta_m6A, mean_mod_delta_psi))
vec_delta_ir, vec_delta_m6A, vec_delta_psi = np.vstack(delta_ir_avg_mod_deltas).T

vec_delta_mods = {
    'm6A': vec_delta_m6A,
    'psi': vec_delta_psi
}

plt.figure(figsize=(10, 5))
for mod_ind, this_mod in enumerate(mods):
    plt.subplot(1, 2, mod_ind+1)
    # comp_vals = []
    for comparator in [lambda x: np.less(x, -0.05), lambda x: np.greater_equal(x, 0.05)]:
        mask = comparator(vec_delta_ir)
        vals = [x for x in vec_delta_mods[this_mod][mask] if ~np.isnan(x)]
        plt.hist(vals, range=[-30, 30], bins=30, alpha=0.5)
    # plt.violinplot(comp_vals, [-1, 1])


# plt.figure(figsize=(10, 5))
# plt.subplot(1, 2, 1)
# plt.plot(vec_delta_ir, vec_delta_m6A, '.')
# plt.xlim([-0.25, 0.25])
# plt.ylim([-30, 30])
# plt.subplot(1, 2, 2)
# plt.plot(vec_delta_ir, vec_delta_psi, '.')
# plt.xlim([-0.25, 0.25])
# plt.ylim([-30, 30])

# delta_ir_max = 0.1
# num_bins = 2
#
# bin_edges = np.linspace(-delta_ir_max, delta_ir_max, num_bins+1)
# bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
# binned_delta_m6A = []
# binned_delta_psi = []
# for bin_i in range(len(bin_centers)):
#     bin_start = bin_edges[bin_i]
#     bin_end = bin_edges[bin_i+1]
#     mask = (vec_delta_ir >= bin_start) * (vec_delta_ir < bin_end)
#     binned_delta_m6A.append(np.nanmean(vec_delta_m6A[mask]))
#     binned_delta_psi.append(np.nanmean(vec_delta_psi[mask]))
#
# binned_deltas = {
#     'm6A': binned_delta_m6A,
#     'psi': binned_delta_psi
# }
#
# xlim = [-delta_ir_max, delta_ir_max]
# xticks = np.linspace(*xlim, num_bins+1)
# plt.figure(figsize=(10, 5))
# for mod_ind, this_mod in enumerate(mods):
#     plt.subplot(1, 2, mod_ind+1)
#     plt.plot(bin_centers, binned_deltas[this_mod], '-o')
#     plt.xlim(xlim)
#     plt.xticks(xticks)
#     plt.xlabel(r'$\Delta$ IR ratio')
#     plt.ylabel(r'Mean $\Delta S$')
#     plt.title(rf'${{{dict_mod_display[this_mod]}}}$')
#     # plt.ylim([-10, 10])
# plt.suptitle(f'{cond_names[1]} vs {cond_names[0]}')