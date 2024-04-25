import os
from glob import glob
import pandas as pd
from gppy import gtf
import numpy as np
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

ds_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT'
gtf_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf'
annot = gtf.parse_gtf(gtf_file)
img_out = 'home/adrian/img_out/transcript_profile'
os.makedirs(img_out, exist_ok=True)
thresh_conf = 80.0
bin_splits = [10, 80, 10]
num_bins = sum(bin_splits)

def get_binned_mod_ratios(in_mod, in_tx):
    start = in_tx.tx_start
    end = in_tx.tx_end
    strand = in_tx.gene.strand

    exon_ranges = [(interval.start, interval.end) for interval in in_tx.exons]
    cds_ranges = [(interval.start, interval.end) for interval in in_tx.cdss]
    five_prime_utr_ranges = [(interval.start, interval.end) for interval in in_tx.five_prime_utrs]
    three_prime_utr_ranges = [(interval.start, interval.end) for interval in in_tx.three_prime_utrs]

    if strand=='-':
        exon_ranges = [this_range[::-1] for this_range in exon_ranges[::-1]]
        cds_ranges = [this_range[::-1] for this_range in cds_ranges[::-1]]
        five_prime_utr_ranges = [this_range[::-1] for this_range in five_prime_utr_ranges[::-1]]
        three_prime_utr_ranges = [this_range[::-1] for this_range in three_prime_utr_ranges[::-1]]

    # print(five_prime_utr_ranges, cds_ranges, three_prime_utr_ranges)

    five_prime_len = sum([max(pair)-min(pair)+ 1 for pair in five_prime_utr_ranges])
    cds_len = sum([max(pair)-min(pair)+ 1 for pair in cds_ranges])
    three_prime_len = sum([max(pair)-min(pair)+ 1 for pair in three_prime_utr_ranges])

    dividers = [1, five_prime_len, five_prime_len + cds_len, five_prime_len + cds_len + three_prime_len]
    bins = [np.linspace(dividers[i], dividers[i+1], this_num_bins+1) for i, this_num_bins in enumerate(bin_splits)]
    bin_edges = np.round(np.concatenate([bins[0][:-1], bins[1][:-1], bins[2]]))

    mod_filtered = in_mod[
        (in_mod['confidence'] >= thresh_conf)
        * (in_mod['chromEnd'] >= min(start, end))
        * (in_mod['chromEnd'] <= max(start, end))
    ]

    binned_mod_ratios = {}
    for mod in ['m6A', 'psi']:
        sub_df = mod_filtered[mod_filtered['name']==mod]
        tx_pos = np.array([in_tx.gpos_to_tpos(gpos)[0] for gpos in sub_df['chromEnd']])
        mod_ratios = sub_df['modRatio'].values
        binned_mod_ratios[mod] = [[] for i in range(len(bin_edges)-1)]
        for i in range(len(bin_edges)-1):
            start_bin = bin_edges[i]
            end_bin = bin_edges[i+1]
            mask = (tx_pos>=start_bin) * (tx_pos<end_bin)
            binned_mod_ratios[mod][i] = list(mod_ratios[mask])

    return binned_mod_ratios, exon_ranges, cds_ranges, five_prime_utr_ranges, three_prime_utr_ranges

def plot_transcript_profile(in_tx, in_mod, exon_ranges, cds_ranges, five_prime_utr_ranges, three_prime_utr_ranges):
    mod_colors = {
        'm6A': 'r',
        'psi': 'g'
    }

    gene_name = in_tx.gene.gene_name
    tx_id = in_tx.tx_id
    start = in_tx.tx_start
    end = in_tx.tx_end

    fig = plt.figure(figsize=(10, 4))

    for mod in ['m6A', 'psi']:
        vec_pos, vec_modRatio = in_mod[in_mod['name']==mod][['chromEnd', f'modRatio']].values.T
        plt.plot(vec_pos, vec_modRatio, f'{mod_colors[mod]}.', label=mod)

    for exon in exon_ranges:
        plt.axvspan(xmin=exon[0], xmax=exon[1], color='grey', alpha=0.1)
    for cds in cds_ranges:
        plt.axvspan(xmin=cds[0], xmax=cds[1], color='blue', alpha=0.1)
    for five_prime_utr in five_prime_utr_ranges:
        plt.axvspan(xmin=five_prime_utr[0], xmax=five_prime_utr[1], color='yellow', alpha=0.1)
    for three_prime_utr in three_prime_utr_ranges:
        plt.axvspan(xmin=three_prime_utr[0], xmax=three_prime_utr[1], color='purple', alpha=0.1)

    plt.title(f'{gene_name}, {tx_id}', fontsize=15)
    plt.xlim([start, end])
    plt.ylim([-5, 105])
    plt.xlabel('Genomic coord.', fontsize=12)
    plt.ylabel('Mod. Ratio', fontsize=12)
    plt.legend(loc='upper left')

    return fig

tx_mod_profile = {}
for this_chr in [str(i) for i in range(1, 23)] + ['X']:
    print(f'chr{this_chr}')
    tx_bed_files = glob(os.path.join(ds_dir, f'chr{this_chr}/transcript', '*.bed'))
    for this_tx_bed_file in tqdm(tx_bed_files):
        geneId, txId = os.path.basename(this_tx_bed_file).rstrip('.bed').split('_')
        this_tx = annot[txId]
        this_tx_bed = pd.read_csv(this_tx_bed_file, sep='\t')
        this_tx_bed_thresh = this_tx_bed[this_tx_bed['confidence']>=thresh_conf]

        tx_binned_mod_ratios, tx_exon_ranges, tx_cds_ranges, tx_five_prime_utr_ranges, tx_three_prime_utr_ranges = \
            get_binned_mod_ratios(this_tx_bed_thresh, this_tx)

        if (len(tx_five_prime_utr_ranges)>0) and (len(tx_cds_ranges)>0) and (len(tx_three_prime_utr_ranges)>0):
            tx_mod_profile[txId] = tx_binned_mod_ratios

        # this_fig = plot_transcript_profile(
        #     this_tx, this_tx_bed_thresh,
        #     tx_exon_ranges, tx_cds_ranges, tx_five_prime_utr_ranges, tx_three_prime_utr_ranges)

agg_m6A_profile = [[] for i in range(num_bins)]
agg_psi_profile = [[] for i in range(num_bins)]
for tx_id, mod_profiles in tx_mod_profile.items():
    for i in range(num_bins):
        agg_m6A_profile[i].extend(mod_profiles['m6A'][i])
        agg_psi_profile[i].extend(mod_profiles['psi'][i])

avg_m6A_profile = [np.mean(x) for x in agg_m6A_profile]
avg_psi_profile = [np.mean(x) for x in agg_psi_profile]
# avg_m6A_profile = [np.median(x) for x in avg_m6A_profile]
# avg_psi_profile = [np.median(x) for x in avg_psi_profile]
# avg_m6A_profile = [np.mean(np.array(x)>50.0) for x in agg_m6A_profile]
# avg_psi_profile = [np.mean(np.array(x)>50.0) for x in agg_psi_profile]

tick_locations = [5, 50, 95]

plt.figure(figsize=(10, 10))
plt.subplot(2, 1, 1)
plt.plot(avg_m6A_profile)
plt.xticks(tick_locations, ['5\' UTR', 'CDS', '3\' UTR'])
plt.tick_params(bottom = False)
plt.axvspan(10, 90, color='blue', alpha=0.1)
plt.xlim([0, 100])
# plt.ylim([0, 20])
plt.ylabel(r'$\langle S \rangle$', fontsize=12)
plt.title('$m^6A$', fontsize=15)
plt.subplot(2, 1, 2)
plt.plot(avg_psi_profile)
plt.xticks(tick_locations, ['5\' UTR', 'CDS', '3\' UTR'])
plt.tick_params(bottom = False)
plt.axvspan(10, 90, color='blue', alpha=0.1)
plt.xlim([0, 100])
# plt.ylim([0, 10])
plt.xlabel('Transcript region', fontsize=12)
plt.ylabel(r'$\langle S \rangle$', fontsize=12)
plt.title('$\psi$', fontsize=15)
plt.suptitle(f'{len(tx_mod_profile)} transcripts', fontsize=15)
plt.savefig(os.path.join(img_out, 'avg_profile.png'), bbox_inches='tight')
plt.close('all')
