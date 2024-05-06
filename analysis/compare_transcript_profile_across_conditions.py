import os
from glob import glob
import pandas as pd
pd.set_option('display.max_columns', None)
from gppy import gtf
import numpy as np
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

### in-house HEK293 ###
ctrl = 'WT'
ctrl_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT'
cond = 'Mettl3-KO'
cond_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293T-Mettl3-KO/rep1'

### NanoSPA ###
# ctrl = 'siCtrl'
# ctrl_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/NanoSPA/HEK_{ctrl}_input_merged'
# # cond = 'siMETTL3'
# cond = 'siTRUB1'
# cond_dir = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/NanoSPA/HEK_{cond}_input_merged'

gtf_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf'
annot = gtf.parse_gtf(gtf_file)
img_out = f'/home/adrian/img_out/transcript_profile_{cond}_{ctrl}'
os.makedirs(img_out, exist_ok=True)
thresh_conf = 80.0
bin_splits = [10, 80, 10]
num_bins = sum(bin_splits)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def get_exon_cds_utr_ranges(in_tx):
    start = in_tx.tx_start
    end = in_tx.tx_end
    strand = in_tx.gene.strand

    exon_ranges = [(interval.start, interval.end) for interval in in_tx.exons]
    five_prime_utr_ranges = [(interval.start, interval.end) for interval in in_tx.five_prime_utrs]
    cds_ranges = [(interval.start, interval.end) for interval in in_tx.cdss]
    three_prime_utr_ranges = [(interval.start, interval.end) for interval in in_tx.three_prime_utrs]

    if strand=='-':
        exon_ranges = [this_range[::-1] for this_range in exon_ranges[::-1]]
        five_prime_utr_ranges = [this_range[::-1] for this_range in five_prime_utr_ranges[::-1]]
        cds_ranges = [this_range[::-1] for this_range in cds_ranges[::-1]]
        three_prime_utr_ranges = [this_range[::-1] for this_range in three_prime_utr_ranges[::-1]]

    five_prime_len = sum([max(pair)-min(pair)+ 1 for pair in five_prime_utr_ranges])
    cds_len = sum([max(pair)-min(pair)+ 1 for pair in cds_ranges])
    three_prime_len = sum([max(pair)-min(pair)+ 1 for pair in three_prime_utr_ranges])

    return {'tx_start_end': (start, end),
            'exon': exon_ranges,
            'five_prime_utr': five_prime_utr_ranges,
            'cds': cds_ranges,
            'three_prime_utr': three_prime_utr_ranges,
            'five_prime_len': five_prime_len,
            'cds_len': cds_len,
            'three_prime_len': three_prime_len,
            }

def get_binned_variable(in_mod, in_tx, dict_ranges, key):
    # print(five_prime_utr_ranges, cds_ranges, three_prime_utr_ranges)
    start_end = dict_ranges['tx_start_end']
    five_prime_len = dict_ranges['five_prime_len']
    cds_len = dict_ranges['cds_len']
    three_prime_len = dict_ranges['three_prime_len']

    dividers = [1, five_prime_len, five_prime_len + cds_len, five_prime_len + cds_len + three_prime_len]
    bins = [np.linspace(dividers[i], dividers[i+1], this_num_bins+1) for i, this_num_bins in enumerate(bin_splits)]
    bin_edges = np.round(np.concatenate([bins[0][:-1], bins[1][:-1], bins[2]]))

    binned_variable = {}
    for mod in ['m6A', 'psi']:
        sub_df = in_mod[in_mod['name']==mod]
        tx_pos = np.array([in_tx.gpos_to_tpos(gpos)[0] for gpos in sub_df['chromEnd']])
        mod_ratios = sub_df[key].values
        binned_variable[mod] = [[] for i in range(len(bin_edges)-1)]
        for i in range(len(bin_edges)-1):
            start_bin = bin_edges[i]
            end_bin = bin_edges[i+1]
            mask = (tx_pos>=start_bin) * (tx_pos<end_bin)
            binned_variable[mod][i] = list(mod_ratios[mask])

    return binned_variable

def plot_transcript_profile(in_tx, in_mod, dict_ranges, key, ylim=None):
    mod_colors = {
        'm6A': 'r',
        'psi': 'g'
    }

    gene_name = in_tx.gene.gene_name
    tx_id = in_tx.tx_id
    start_end = dict_ranges['tx_start_end']

    fig = plt.figure(figsize=(10, 8))

    plt.subplot(2, 1, 1)
    for mod in ['m6A', 'psi']:
        vec_pos, vec_variable = in_mod[in_mod['name']==mod][['chromEnd', key]].values.T
        plt.plot(vec_pos, vec_variable, f'{mod_colors[mod]}.', label=mod)
    for exon in dict_ranges['exon']:
        plt.axvspan(xmin=exon[0], xmax=exon[1], color='grey', alpha=0.1)
    for cds in dict_ranges['cds']:
        plt.axvspan(xmin=cds[0], xmax=cds[1], color='blue', alpha=0.1)
    for five_prime_utr in dict_ranges['five_prime_utr']:
        plt.axvspan(xmin=five_prime_utr[0], xmax=five_prime_utr[1], color='yellow', alpha=0.1)
    for three_prime_utr in dict_ranges['three_prime_utr']:
        plt.axvspan(xmin=three_prime_utr[0], xmax=three_prime_utr[1], color='purple', alpha=0.1)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlim(start_end)
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel('Genomic coord.', fontsize=12)
    plt.ylabel(f'$S_{{{cond}}} - S_{{{ctrl}}}$', fontsize=12)
    plt.legend(loc='upper left')

    five_prime_len = dict_ranges['five_prime_len']
    cds_len = dict_ranges['cds_len']
    three_prime_len = dict_ranges['three_prime_len']
    tx_len = five_prime_len + cds_len + three_prime_len

    plt.subplot(2, 1, 2)
    for mod in ['m6A', 'psi']:
        vec_pos, vec_variable = in_mod[in_mod['name']==mod][['tpos', key]].values.T
        plt.plot(vec_pos, vec_variable, f'{mod_colors[mod]}-o', label=mod)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.axvspan(five_prime_len, five_prime_len + cds_len, color='blue', alpha=0.1)
    plt.xlim([0, tx_len])
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel('Transcriptomic coord.', fontsize=12)
    plt.ylabel(f'$S_{{{cond}}} - S_{{{ctrl}}}$', fontsize=12)
    tick_locations = [five_prime_len*0.5, five_prime_len+cds_len*0.5, five_prime_len+cds_len+three_prime_len*0.5]
    plt.xticks(tick_locations, ['5\' UTR', 'CDS', '3\' UTR'])
    plt.tick_params(bottom=False)

    plt.suptitle(f'{gene_name}, {tx_id}', fontsize=15)

    return fig


def get_merged_bed(in_ctrl_bed, in_cond_bed, tx_annot, dict_ranges):
    start_end = dict_ranges['tx_start_end']

    merged_bed = pd.merge(in_ctrl_bed, in_cond_bed,
                          on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                          suffixes=[f'_{ctrl}', f'_{cond}'])
    merged_bed = merged_bed[
        (merged_bed['chromEnd'] >= min(start_end))
        * (merged_bed['chromEnd'] <= max(start_end))
        ]
    merged_bed['tpos'] = [tx_annot.gpos_to_tpos(gpos)[0] for gpos in merged_bed['chromEnd']]
    merged_bed['delta'] = merged_bed[f'modRatio_{cond}'] - merged_bed[f'modRatio_{ctrl}']

    return merged_bed


def get_r_region_delta_tuples_per_exon(in_bed, in_dict_ranges, primary_mod, secondary_mod):
    exon_r_region_delta_tuples = []
    for this_exon in in_dict_ranges['exon']:
        exon_bed = in_bed[(in_bed['chromEnd']>=min(this_exon)) * (in_bed['chromEnd']<=max(this_exon))]
        primary_tpos, primary_delta, primary_region = exon_bed[exon_bed['name']==primary_mod][['tpos', 'delta', 'region']].values.T
        secondary_tpos, secondary_delta, secondary_region = exon_bed[exon_bed['name']==secondary_mod][['tpos', 'delta', 'region']].values.T
        r_mat = primary_tpos[:, np.newaxis] - secondary_tpos[np.newaxis, :]
        region_pairs = [(prim_region, sec_region) for prim_region in primary_region for sec_region in secondary_region]
        delta_pairs = [(prim_delta, sec_delta) for prim_delta in primary_delta for sec_delta in secondary_delta]
        exon_r_region_delta_tuples.extend([(r, *region_pairs[i], *delta_pairs[i]) for (i, r) in enumerate(r_mat.flatten())])
    return exon_r_region_delta_tuples


def get_binned_delta_by_dist(in_r_region_delta_tuples, div_pts, bin_edges, prim_mod, sec_mod, prim_region=None, sec_region=None):
    # r_region_delta_tuples = np.vstack([sublist for tx_list in in_r_region_delta_tuples.values() for sublist in tx_list])
    # r_delta_tuples = np.float64(r_region_delta_tuples[
    #                      (r_region_delta_tuples[:, 1]==prim_region)
    #                      * (r_region_delta_tuples[:, 2]==sec_region)][:, [0, 3, 4]])
    if (prim_region is None) and (sec_region is None):
        r_delta_tuples = [(tup[0], tup[3], tup[4]) for tx_list in in_r_region_delta_tuples.values() for tup in tx_list]
    else:
        r_delta_tuples = [(tup[0], tup[3], tup[4]) for tx_list in in_r_region_delta_tuples.values()
                                    for tup in tx_list if ((tup[1]==prim_region) and (tup[2]==sec_region))]
    if len(r_delta_tuples):
        r_delta_tuples = np.vstack(r_delta_tuples)
    else:
        return None, None

    dict_col = {
        'm6A': 1,
        'psi': 2
    }

    stratified_r_binned_delta = []
    stratified_num_samples = []
    for prim_delta_start, prim_delta_end in zip(div_pts[:-1], div_pts[1:]):
        r_delta = r_delta_tuples[(r_delta_tuples[:, dict_col[prim_mod]]>=prim_delta_start) * (r_delta_tuples[:, dict_col[prim_mod]]<prim_delta_end)]
        binned_r_delta = [[] for r in bin_edges[:-1]]
        num_samples = 0
        for r_ind, (r_start, r_end) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
            r_subsample = r_delta[(r_delta[:, 0]>=r_start) * (r_delta[:, 0]<r_end)][:, dict_col[sec_mod]]
            num_samples += len(r_subsample)
            binned_r_delta[r_ind] = np.mean(r_subsample)
            # binned_r_delta[r_ind] = np.median(r_subsample)
        stratified_r_binned_delta.append(binned_r_delta)
        stratified_num_samples.append(num_samples)

    return stratified_r_binned_delta, stratified_num_samples


def annotate_region(in_bed, in_ranges):
    out_bed = in_bed.copy()
    out_bed['region'] = 'None'
    for region in ['five_prime_utr', 'cds', 'three_prime_utr']:
        for annot in in_ranges[region]:
            out_bed.loc[(in_bed['chromEnd'] <= max(annot)) * (in_bed['chromEnd'] >= min(annot)), 'region'] = region
    return out_bed


def plot_stratified_corr(strat_div_pts, strat_r_binned_delta, bin_centers, strat_num_samples, prim_mod, sec_mod, title=None):
    stratified_labels = [f'[{start},{end})' for start, end in zip(strat_div_pts[:-1], strat_div_pts[1:])]
    plt.figure(figsize=(16, 9))
    for strat, label, n_samples in zip(strat_r_binned_delta, stratified_labels, strat_num_samples):
        plt.plot(strat, label=f'{label}, {n_samples}')
        plt.xticks(range(len(bin_centers)), bin_centers)
    plt.legend(title=f'$\Delta S_{dict_mod_display[prim_mod]}$', loc='lower left')
    plt.xlabel(f'Distance from ${dict_mod_display[prim_mod]}$ site', fontsize=12)
    plt.ylabel(f'$\Delta S_{dict_mod_display[sec_mod]}$', fontsize=12)
    plt.ylim([-15, 15])
    plt.axhline(y=0, color='gray', linestyle='--')
    if title is not None:
        plt.title(title, fontsize=15)
    plt.show()

# tx_diff_profile = {}
tx_merged_bed = {}
tx_r_region_delta_tuples = {}
tx_deltas = {}
for this_chr in [str(i) for i in range(1, 23)] + ['X']:
    print(f'chr{this_chr}')
    ctrl_bed_files = glob(os.path.join(ctrl_dir, f'chr{this_chr}/transcript', '*.bed'))
    cond_bed_files = glob(os.path.join(cond_dir, f'chr{this_chr}/transcript', '*.bed'))
    ctrl_geneId_txId = [os.path.basename(f).rstrip('.bed') for f in ctrl_bed_files]
    cond_geneId_txId = [os.path.basename(f).rstrip('.bed') for f in cond_bed_files]
    common_geneId_txId = list(set(ctrl_geneId_txId).intersection(set(cond_geneId_txId)))
    common_geneId_txId.sort()

    for this_common_geneId_txId in tqdm(common_geneId_txId):
        geneId, txId = this_common_geneId_txId.split('_')
        this_tx_annot = annot[txId]

        ctrl_bed = pd.read_csv(os.path.join(ctrl_dir, f'chr{this_chr}/transcript', f'{this_common_geneId_txId}.bed'), sep='\t')
        cond_bed = pd.read_csv(os.path.join(cond_dir, f'chr{this_chr}/transcript', f'{this_common_geneId_txId}.bed'), sep='\t')

        ctrl_bed_thresh = ctrl_bed[ctrl_bed['confidence']>=thresh_conf]
        cond_bed_thresh = cond_bed[cond_bed['confidence']>=thresh_conf]

        this_tx_dict_ranges = get_exon_cds_utr_ranges(this_tx_annot)
        if ((len(this_tx_dict_ranges['five_prime_utr'])==0)
                or (len(this_tx_dict_ranges['cds'])==0)
                or (len(this_tx_dict_ranges['three_prime_utr'])==0)):
            continue

        this_tx_merged_bed = get_merged_bed(ctrl_bed_thresh, cond_bed_thresh, this_tx_annot, this_tx_dict_ranges)
        this_tx_merged_bed = annotate_region(this_tx_merged_bed, this_tx_dict_ranges)
        tx_merged_bed[txId] = this_tx_merged_bed

        tx_deltas[txId] = {}
        tx_deltas[txId]['m6A'] = this_tx_merged_bed[this_tx_merged_bed['name'] == 'm6A']['delta'].values
        tx_deltas[txId]['psi'] = this_tx_merged_bed[this_tx_merged_bed['name'] == 'psi']['delta'].values

        # this_fig = plot_transcript_profile(this_tx_annot, this_tx_merged_bed, this_tx_dict_ranges, key='delta')
        # this_fig.savefig(os.path.join(img_out, f'{this_common_geneId_txId}.png'), bbox_inches='tight')
        # plt.close('all')

        # tx_binned_deltas = get_binned_variable(this_tx_merged_bed, this_tx_annot, this_tx_start_end, this_tx_dict_ranges, key='delta')
        # tx_diff_profile[txId] = tx_binned_deltas

        tx_r_region_delta_tuples[txId] = get_r_region_delta_tuples_per_exon(this_tx_merged_bed, this_tx_dict_ranges, 'm6A', 'psi')

r_bin_width = 10
r_bin_max = 105
r_bin_edges = np.arange(-r_bin_max, r_bin_max + r_bin_width, r_bin_width)
r_bin_centers = (r_bin_edges[1:] + r_bin_edges[:-1]) // 2

primary_div_pts = np.array([-50, -30, -10, 10, 30])

binned_delta, num_samples = get_binned_delta_by_dist(tx_r_region_delta_tuples, primary_div_pts, r_bin_edges, prim_mod='m6A', sec_mod='psi')
plot_stratified_corr(primary_div_pts, binned_delta, r_bin_centers, num_samples, prim_mod='m6A', sec_mod='psi')

binned_delta, num_samples = get_binned_delta_by_dist(tx_r_region_delta_tuples, primary_div_pts, r_bin_edges, prim_mod='psi', sec_mod='m6A')
plot_stratified_corr(primary_div_pts, binned_delta, r_bin_centers, num_samples, prim_mod='psi', sec_mod='m6A')

# for region in ['cds', 'three_prime_utr']:
#     binned_delta, num_samples = get_binned_delta_by_dist(
#         tx_r_region_delta_tuples, primary_div_pts, r_bin_edges, region, region
#     )
#     plot_stratified_corr(primary_div_pts, binned_delta, r_bin_edges, num_samples, title=region)

########################################################################################################################
# all_m6A_deltas = np.concatenate([this_delta['m6A'] for this_delta in tx_deltas.values()])
# all_psi_deltas = np.concatenate([this_delta['psi'] for this_delta in tx_deltas.values()])
#
# plt.figure(figsize=(5, 5))
# plt.hist(all_m6A_deltas, bins=40, range=[-20, 20])
# plt.xlabel('$\Delta S_{m^6A}$', fontsize=12)
# plt.ylabel('Num. sites', fontsize=12)
# # plt.yscale('log')
#
# plt.figure(figsize=(5, 5))
# plt.hist(all_psi_deltas, bins=40, range=[-20, 20])
# plt.xlabel('$\Delta S_{\psi}$', fontsize=12)
# plt.ylabel('Num. sites', fontsize=12)
# # plt.yscale('log')

tx_avg_deltas = {
    tx: (np.mean(deltas['m6A']), np.mean(deltas['psi'])) for tx, deltas in tx_deltas.items()
    if (len(deltas['m6A'])>=10 and len(deltas['psi'])>=10)
}
tx_avg_delta_m6A, tx_avg_delta_psi = np.vstack(list(tx_avg_deltas.values())).T

# plt.figure(figsize=(5, 5))
# plt.hist(tx_avg_delta_m6A, bins=35, range=[-30, 5])
# plt.xlabel(r'$\langle \Delta S_{m^6A} \rangle$ per transcript', fontsize=12)
# plt.ylabel('Num. transcripts', fontsize=12)
# plt.axvline(x=0, c='r', linestyle='--')
#
# plt.figure(figsize=(5, 5))
# plt.hist(tx_avg_delta_psi, bins=12, range=[-6, 6])
# plt.xlabel(r'$\langle \Delta S_{\psi} \rangle$ per transcript', fontsize=12)
# plt.ylabel('Num. transcripts', fontsize=12)
# plt.axvline(x=0, c='r', linestyle='--')

xlim = [-30, 5]
ylim = [-6, 6]
xbins = xlim[1] - xlim[0]
ybins = (ylim[1] - ylim[0]) * 2

fig = plt.figure(figsize=(10, 10), layout='constrained')
gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

ax.plot(tx_avg_delta_m6A, tx_avg_delta_psi, '.')
ax.axvline(x=0, c='r', linestyle='--')
ax.axhline(y=0, c='r', linestyle='--')
ax.set_xlabel(r'$\langle \Delta S_{m^6A} \rangle$ per transcript', fontsize=12)
ax.set_ylabel(r'$\langle \Delta S_{\psi} \rangle$ per transcript', fontsize=12)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax_histx.hist(tx_avg_delta_m6A, bins=xbins, range=xlim)
ax_histx.axvline(x=0, c='r', linestyle='--')
ax_histx.set_xlim(xlim)
# ax_histx.set_ylim(ylim)
ax_histx.set_ylabel('Num. transcripts', fontsize=12)

ax_histy.hist(tx_avg_delta_psi, bins=ybins, range=ylim, orientation='horizontal')
ax_histy.axhline(y=0, c='r', linestyle='--')
# ax_histy.set_xlim(xlim)
ax_histy.set_ylim(ylim)
ax_histy.set_xlabel('Num. transcripts', fontsize=12)

fig.savefig(os.path.join(img_out, 'avg_delta_per_transcript.png'), bbox_inches='tight')

########################################################################################################################
min_sites = 10

ctrl_tx_avg_mod_ratio = {
    this_tx: (this_bed[this_bed['name']=='m6A'][f'modRatio_{ctrl}'].mean(),
              this_bed[this_bed['name']=='psi'][f'modRatio_{ctrl}'].mean())
    for this_tx, this_bed in tx_merged_bed.items()
    if (len(this_bed[this_bed['name'] == 'm6A']) >= min_sites) and (len(this_bed[this_bed['name'] == 'm6A']) >= min_sites)
}
ctrl_tx_avg_mod_ratio_m6A, ctrl_tx_avg_mod_ratio_psi = np.vstack(list(ctrl_tx_avg_mod_ratio.values())).T
ctrl_tx_avg_mod_ratio_m6A = np.nan_to_num(ctrl_tx_avg_mod_ratio_m6A)
ctrl_tx_avg_mod_ratio_psi = np.nan_to_num(ctrl_tx_avg_mod_ratio_psi)

cond_tx_avg_mod_ratio = {
    this_tx: (this_bed[this_bed['name']=='m6A'][f'modRatio_{cond}'].mean(),
              this_bed[this_bed['name']=='psi'][f'modRatio_{cond}'].mean())
    for this_tx, this_bed in tx_merged_bed.items()
    if (len(this_bed[this_bed['name'] == 'm6A']) >= min_sites) and (len(this_bed[this_bed['name'] == 'm6A']) >= min_sites)
}
cond_tx_avg_mod_ratio_m6A, cond_tx_avg_mod_ratio_psi = np.vstack(list(cond_tx_avg_mod_ratio.values())).T
cond_tx_avg_mod_ratio_m6A = np.nan_to_num(cond_tx_avg_mod_ratio_m6A)
cond_tx_avg_mod_ratio_psi = np.nan_to_num(cond_tx_avg_mod_ratio_psi)


xlim = [-5, 50]
ylim = [-1, 25]

plt.figure(figsize=(10, 10))

plt.subplot(2, 2, 1)
plt.plot(ctrl_tx_avg_mod_ratio_m6A, ctrl_tx_avg_mod_ratio_psi, '.')
plt.xlim(xlim)
plt.ylim(ylim)
# plt.xlabel(r'$\langle S_{m^6A} \rangle$ per transcript', fontsize=12)
plt.ylabel(r'$\langle S_{\psi} \rangle$ per transcript', fontsize=12)
plt.title(ctrl, fontsize=15)

plt.subplot(2, 2, 2)
plt.plot(cond_tx_avg_mod_ratio_m6A, cond_tx_avg_mod_ratio_psi, '.')
plt.xlim(xlim)
plt.ylim(ylim)
# plt.xlabel(r'$\langle S_{m^6A} \rangle$ per transcript', fontsize=12)
# plt.ylabel(r'$\langle S_{\psi} \rangle$ per transcript', fontsize=12)
plt.title(cond, fontsize=15)

tracked_ind = {
    'm6A': 0,
    'psi': 1
}
tracked_mod = 'm6A'
tracked_mod_thresh = 15.0

above_track_tx = [tx for tx, avg_mod_ratio in ctrl_tx_avg_mod_ratio.items() if avg_mod_ratio[tracked_ind[tracked_mod]]>=tracked_mod_thresh]
plt.subplot(2, 2, 3)
for ind, tx in enumerate(above_track_tx):
    if ind==0:
        plt.plot(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], 'b.', label=ctrl)
        plt.plot(cond_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][1], 'r.', label=cond)
    else:
        plt.plot(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], 'b.')
        plt.plot(cond_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][1], 'r.')
    plt.plot((ctrl_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][0]),
             (ctrl_tx_avg_mod_ratio[tx][1], cond_tx_avg_mod_ratio[tx][1]), c='grey', linestyle='-', alpha=0.5)
    # dx = cond_tx_avg_mod_ratio[tx][0] - ctrl_tx_avg_mod_ratio[tx][0]
    # dy = cond_tx_avg_mod_ratio[tx][1] - ctrl_tx_avg_mod_ratio[tx][1]
    # plt.arrow(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], dx, dy,
    #           color='grey', alpha=0.5, head_width=1)
plt.legend(loc='upper right')
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel(r'$\langle S_{m^6A} \rangle$ per transcript', fontsize=12)
plt.ylabel(r'$\langle S_{\psi} \rangle$ per transcript', fontsize=12)
plt.title(fr'$\langle S_{{{dict_mod_display[tracked_mod]}}}^{{{ctrl}}} \rangle \geq {tracked_mod_thresh}$', fontsize=15)

tracked_mod = 'psi'
tracked_mod_thresh = 10.0

above_track_tx = [tx for tx, avg_mod_ratio in ctrl_tx_avg_mod_ratio.items() if avg_mod_ratio[tracked_ind[tracked_mod]]>=tracked_mod_thresh]
plt.subplot(2, 2, 4)
for ind, tx in enumerate(above_track_tx):
    if ind==0:
        plt.plot(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], 'b.', label=ctrl)
        plt.plot(cond_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][1], 'r.', label=cond)
    else:
        plt.plot(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], 'b.')
        plt.plot(cond_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][1], 'r.')
    plt.plot((ctrl_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][0]),
             (ctrl_tx_avg_mod_ratio[tx][1], cond_tx_avg_mod_ratio[tx][1]), c='grey', linestyle='-', alpha=0.5)
    # dx = cond_tx_avg_mod_ratio[tx][0] - ctrl_tx_avg_mod_ratio[tx][0]
    # dy = cond_tx_avg_mod_ratio[tx][1] - ctrl_tx_avg_mod_ratio[tx][1]
    # plt.arrow(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], dx, dy,
    #           color='grey', alpha=0.5, head_width=1)
plt.legend(loc='upper right')
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel(r'$\langle S_{m^6A} \rangle$ per transcript', fontsize=12)
# plt.ylabel(r'$\langle S_{\psi} \rangle$ per transcript', fontsize=12)
plt.title(fr'$\langle S_{{{dict_mod_display[tracked_mod]}}}^{{{ctrl}}} \rangle \geq {tracked_mod_thresh}$', fontsize=15)

plt.savefig(os.path.join(img_out, f'avg_S_per_transcript_{ctrl}_{cond}.png'), bbox_inches='tight')

# below_track_tx = [tx for tx, avg_mod_ratio in ctrl_tx_avg_mod_ratio.items() if avg_mod_ratio[0]<tracked_mod_thresh]
# plt.subplot(2, 2, 4)
# for ind, tx in enumerate(below_track_tx):
#     if ind==0:
#         plt.plot(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], 'b.', label='WT')
#         plt.plot(cond_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][1], 'r.', label='Mettl3-KO')
#     else:
#         plt.plot(ctrl_tx_avg_mod_ratio[tx][0], ctrl_tx_avg_mod_ratio[tx][1], 'b.')
#         plt.plot(cond_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][1], 'r.')
#     plt.plot((ctrl_tx_avg_mod_ratio[tx][0], cond_tx_avg_mod_ratio[tx][0]),
#              (ctrl_tx_avg_mod_ratio[tx][1], cond_tx_avg_mod_ratio[tx][1]), c='grey', linestyle='-', alpha=0.5)
# plt.legend(loc='upper right')
# plt.xlim(xlim)
# plt.ylim(ylim)
# plt.xlabel(r'$\langle S_{m^6A} \rangle$ per transcript', fontsize=12)
# plt.ylabel(r'$\langle S_{\psi} \rangle$ per transcript', fontsize=12)
# plt.title(fr'$\langle S_{{m^6A}} \rangle < {tracked_mod_thresh}$', fontsize=15)

# common_tx = list(set(ctrl_tx_avg_mod_ratio.keys()).intersection(set(cond_tx_avg_mod_ratio.keys())))
# delta_avg_S = {tx: (cond_tx_avg_mod_ratio[tx][0]-ctrl_tx_avg_mod_ratio[tx][0],
#                     cond_tx_avg_mod_ratio[tx][1]-ctrl_tx_avg_mod_ratio[tx][1])
#                for tx in ctrl_tx_avg_mod_ratio.keys()}
# delta_avg_S_m6A, delta_avg_S_psi = np.vstack(list(delta_avg_S.values())).T
# plt.figure(figsize=(5, 5))
# plt.plot(ctrl_tx_avg_mod_ratio_m6A, delta_avg_S_psi, '.')


########################################################################################################################
# import shutil
#
# margin = 2
#
# low_delta_psi_tx = [k for k, v in tx_avg_deltas.items() if v[1]<-margin]
# high_delta_psi_tx = [k for k, v in tx_avg_deltas.items() if v[1]>=margin]
#
# sub_img_out = os.path.join(img_out, 'high_delta_psi_tx')
# os.makedirs(sub_img_out, exist_ok=True)
# for this_tx in high_delta_psi_tx:
#     filename = os.path.basename(glob(os.path.join(img_out, f'*_{this_tx}.png'))[0])
#     shutil.copyfile(os.path.join(img_out, filename), os.path.join(sub_img_out, filename))
#
# sub_img_out = os.path.join(img_out, 'low_delta_psi_tx')
# os.makedirs(sub_img_out, exist_ok=True)
# for this_tx in low_delta_psi_tx:
#     filename = os.path.basename(glob(os.path.join(img_out, f'*_{this_tx}.png'))[0])
#     shutil.copyfile(os.path.join(img_out, filename), os.path.join(sub_img_out, filename))
#
# low_delta_psi_genes = [annot[k].gene.gene_name for k, v in tx_avg_deltas.items() if v[1]<-margin]
# high_delta_psi_genes = [annot[k].gene.gene_name for k, v in tx_avg_deltas.items() if v[1]>=margin]
# common_genes = list(set(low_delta_psi_genes).intersection(set(high_delta_psi_genes)))

# 'SERBP1'
# low_delta_psi_tx = ['ENST00000370994']
# high_delta_psi_tx = ['ENST00000361219', 'ENST00000370995']


# plt.figure(figsize=(5, 5))
# plt.hist(tx_avg_psi_deltas.values(), bins=40, range=[-20, 20])
# plt.xlabel('tx avg. $\Delta S_{\psi}$', fontsize=12)
# plt.ylabel('Num. sites', fontsize=12)

# tx_high_delta_psi = [(tx, avg_delta) for tx, avg_delta in tx_avg_psi_deltas.items() if avg_delta>=5]

########################################################################################################################

# agg_m6A_profile = [[] for i in range(num_bins)]
# agg_psi_profile = [[] for i in range(num_bins)]
# for tx_id, mod_profiles in tx_diff_profile.items():
#     for i in range(num_bins):
#         agg_m6A_profile[i].extend(mod_profiles['m6A'][i])
#         agg_psi_profile[i].extend(mod_profiles['psi'][i])
#
# avg_m6A_profile = [np.mean(x) for x in agg_m6A_profile]
# avg_psi_profile = [np.mean(x) for x in agg_psi_profile]
# # avg_m6A_profile = [np.median(x) for x in avg_m6A_profile]
# # avg_psi_profile = [np.median(x) for x in avg_psi_profile]
# # avg_m6A_profile = [np.mean(np.array(x)>50.0) for x in agg_m6A_profile]
# # avg_psi_profile = [np.mean(np.array(x)>50.0) for x in agg_psi_profile]
#
# tick_locations = [5, 50, 95]
#
# plt.figure(figsize=(10, 10))
# plt.subplot(2, 1, 1)
# plt.plot(avg_m6A_profile)
# plt.xticks(tick_locations, ['5\' UTR', 'CDS', '3\' UTR'])
# plt.tick_params(bottom = False)
# plt.axvspan(10, 90, color='blue', alpha=0.1)
# plt.xlim([0, 100])
# plt.ylim([-10, 5])
# plt.axhline(y=0, color='r', linestyle='--')
# plt.ylabel(rf'$ \langle S_{{{cond}}} - S_{{{ctrl}}} \rangle $', fontsize=12)
# plt.title('$m^6A$', fontsize=15)
# plt.subplot(2, 1, 2)
# plt.plot(avg_psi_profile)
# plt.xticks(tick_locations, ['5\' UTR', 'CDS', '3\' UTR'])
# plt.tick_params(bottom=False)
# plt.axvspan(10, 90, color='blue', alpha=0.1)
# plt.xlim([0, 100])
# plt.ylim([-2, 2])
# plt.axhline(y=0, color='r', linestyle='--')
# plt.xlabel('Transcript region', fontsize=12)
# plt.ylabel(rf'$ \langle S_{{{cond}}} - S_{{{ctrl}}} \rangle $', fontsize=12)
# plt.title('$\psi$', fontsize=15)
# plt.suptitle(f'{len(tx_diff_profile)} transcripts', fontsize=15)
# plt.savefig(os.path.join(img_out, 'avg_profile.png'), bbox_inches='tight')
# plt.close('all')
