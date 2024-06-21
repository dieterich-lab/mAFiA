import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

width = 0.4

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
gene_name = 'Synpo2l'
events = ['3 -> 5', 'AFE']
colors = ['red', 'green']
counts = {
    'TAC': {
        'day1': [198, 11],
        'day7': [166, 135],
        'day21': [227, 468],
        'day56': [197, 471]
    },
    'SHAM': {
        'day1': [161, 10],
        'day7': [119, 19],
        'day21': [170, 10],
        'day56': [117, 47]
    }
}

event_colors = {ev: c for ev, c in zip(events, colors)}

img_out = '/home/adrian/img_out/TAC_splice_site_mod_changes'
os.makedirs(img_out, exist_ok=True)

days = [
    'day1',
    'day7',
    'day21',
    'day56'
]


def normalize_and_recast(in_counts, in_events):
    norm_counts = np.array([[this_count[0]/sum(this_count), this_count[1]/sum(this_count)] for this_count in in_counts.values()])
    return {in_events[i]: norm_counts[:, i] for i in range(len(events))}


fig = plt.figure(figsize=(8, 8))
for subplot_ind, ds in enumerate(['SHAM', 'TAC']):
    ds_ratios = normalize_and_recast(counts[ds], events)
    ax = fig.add_subplot(2, 1, subplot_ind+1)
    bottom = np.zeros(len(days))
    for this_event, this_ratio in ds_ratios.items():
        p = ax.bar(days, this_ratio, width, label=this_event, bottom=bottom, color=event_colors[this_event])
        bottom += this_ratio
    plt.ylabel('Splice Variant Ratio (%)', fontsize=12)
    plt.ylim([0, 1.1])
    plt.yticks(np.arange(0, 1.1, 0.25), np.int32(np.arange(0, 1.1, 0.25)*100))
    plt.legend(loc='upper right', title='Splice variant')
    plt.title(ds, fontsize=12)
plt.xlabel('Days', fontsize=12)
plt.suptitle(gene_name, fontsize=15)
fig.savefig(os.path.join(img_out, f'{gene_name}_AS_progression.png'), bbox_inches='tight')
