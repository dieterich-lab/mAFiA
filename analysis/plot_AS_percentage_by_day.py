import os
import numpy as np
import pysam
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from functools import reduce


bar_width = 0.4


def get_long_read_counts_from_bam(in_res_dir, in_gene_name, in_events):
    out_dict = {}
    for this_cond in ['TAC', 'SHAM']:
        out_dict[this_cond] = {}
        for this_day in ['day1', 'day7', 'day21', 'day56']:
            out_dict[this_cond][this_day] = []
            for this_event in in_events:
                bam_file = os.path.join(in_res_dir, '_'.join([this_cond, this_day]), f'{in_gene_name}.{this_event}.mAFiA.reads.bam')
                out_dict[this_cond][this_day].append(
                    reduce(lambda x, y: x + y,
                           [int(this_line.split('\t')[2]) for this_line in pysam.idxstats(bam_file).split('\n') if
                            len(this_line)])
                )
    return out_dict


def normalize_and_recast(in_counts, in_events):
    norm_counts = np.array([[this_count[i]/sum(this_count) for i in range(len(events))] for this_count in in_counts.values()])
    return {in_events[i]: norm_counts[:, i] for i in range(len(events))}


def plot_bar_chart(in_counts, long_short):
    fig = plt.figure(figsize=(8, 8))
    for subplot_ind, ds in enumerate(['SHAM', 'TAC']):
        ds_ratios = normalize_and_recast(in_counts[ds], events)
        ax = fig.add_subplot(2, 1, subplot_ind+1)
        bottom = np.zeros(len(days))
        for this_event, this_ratio in ds_ratios.items():
            p = ax.bar(days, this_ratio, bar_width, label=dict_event_isoform[this_event], bottom=bottom, color=event_colors[this_event])
            bottom += this_ratio
        total_counts = [sum(this_val) for this_val in in_counts[ds].values()]
        for x_loc, this_total_count in enumerate(total_counts):
            plt.text(x_loc, 1.02, f'({this_total_count})', fontsize=10, horizontalalignment='center')
        plt.ylabel('Splice Variant Ratio (%)', fontsize=12)
        plt.ylim([0, 1.1])
        plt.yticks(np.arange(0, 1.1, 0.25), np.int32(np.arange(0, 1.1, 0.25)*100))
        if subplot_ind==0:
            # plt.legend(loc='lower right', title='Splice variant')
            plt.legend(loc='lower right')
        plt.title(ds, fontsize=12)
    plt.xlabel('Days', fontsize=12)

    if long_short=='long':
        plt.suptitle(f'{gene_name}\nONT', fontsize=15)
        fig.savefig(os.path.join(img_out, f'{gene_name}_AS_progression_long_reads.png'), bbox_inches='tight')
    else:
        plt.suptitle(f'{gene_name}\nIllumina', fontsize=15)
        fig.savefig(os.path.join(img_out, f'{gene_name}_AS_progression_short_reads.png'), bbox_inches='tight')


### Fhl1 ###
# gene_name = 'Fhl1'
# events = ['3 -> 9', '5 -> 9']
# colors = ['red', 'purple']
# counts = {
#     'TAC': {
#         'day1': [119, 25],
#         'day7': [158, 593],
#         'day21': [168, 820],
#         'day56': [189, 865]
#     },
#     'SHAM': {
#         'day1': [109, 14],
#         'day7': [147, 15],
#         'day21': [172, 19],
#         'day56': [153, 63]
#     }
# }

### Synpo2l ###
# gene_name = 'Synpo2l'
# events = ['3 -> 5', 'AFE']
# colors = ['red', 'green']
# counts = {
#     'TAC': {
#         'day1': [198, 11],
#         'day7': [166, 135],
#         'day21': [227, 468],
#         'day56': [197, 471]
#     },
#     'SHAM': {
#         'day1': [161, 10],
#         'day7': [119, 19],
#         'day21': [170, 10],
#         'day56': [117, 47]
#     }
# }

### Rcan1 ###
gene_name = 'Rcan1'
events = ['2 -> 6', '4 -> 6']
isoforms = ['ENSMUST00000060005', 'ENSMUST00000023672']
dict_event_isoform = {ev: iso for (ev, iso) in zip(events, isoforms)}
colors = ['blue', 'purple']
short_counts = {
    'TAC': {
        'day1': [146, 14],
        'day7': [102, 372],
        'day21': [149, 307],
        'day56': [149, 327]
    },
    'SHAM': {
        'day1': [100, 13],
        'day7': [88, 27],
        'day21': [108, 66],
        'day56': [126, 91]
    }
}

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'
long_counts = get_long_read_counts_from_bam(res_dir, gene_name, ['exon2', 'exon4'])


### Vdac3 ###
# gene_name = 'Vdac3'
# events = ['lsv1', 'lsv2', 'lsv3']
# colors = ['red', 'green', 'blue']
# counts = {
#     'TAC': {
#         'day1': [1204, 205, 56],
#         'day7': [1079, 124, 70],
#         'day21': [995, 184, 75],
#         'day56': [737, 1125, 66]
#     },
#     'SHAM': {
#         'day1': [1437, 91, 44],
#         'day7': [1211, 84, 52],
#         'day21': [1327, 151, 64],
#         'day56': [1142, 0, 45]
#     }
# }

### Csad ###
# gene_name = 'Csad'
# events = ['lsv1', 'lsv2', 'lsv3']
# colors = ['blue', 'red', 'green']
# counts = {
#     'TAC': {
#         'day1': [21, 9, 8],
#         'day7': [30, 4, 1],
#         'day21': [25, 5, 1],
#         'day56': [18, 4, 3]
#     },
#     'SHAM': {
#         'day1': [16, 6, 4],
#         'day7': [14, 7, 9],
#         'day21': [15, 10, 4],
#         'day56': [17, 13, 11]
#     }
# }


event_colors = {ev: c for ev, c in zip(events, colors)}

img_out = '/home/adrian/img_out/TAC_splice_site_mod_changes'
os.makedirs(img_out, exist_ok=True)

days = [
    'day1',
    'day7',
    'day21',
    'day56'
]

plot_bar_chart(short_counts, 'short')
plot_bar_chart(long_counts, 'long')