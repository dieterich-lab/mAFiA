import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import os
import numpy as np

def plot_bar_chart_site_trend_by_region(img_out, fname, dict_mod_display):
    regions = ['five_prime_utr', 'CDS', 'three_prime_utr']
    fig = plt.figure(figsize=(10, 5))
    all_region_site_counts = []
    for subplot_ind, mod in enumerate(['psi', 'm6A']):
        plt.subplot(1, 2, subplot_ind + 1)
        for mask_name in ['increasing', 'decreasing']:
            this_df = pd.read_csv(
                os.path.join(img_out, fname.format(mask_name)),
                sep='\t')
            all_counts = Counter(this_df[this_df['name'] == mod]['region'])
            region_site_counts = [all_counts[this_region] for this_region in regions]
            all_region_site_counts.extend(region_site_counts)
            if mask_name == 'decreasing':
                region_site_counts = [-this_count for this_count in region_site_counts]
                color = 'cyan'
            else:
                color = 'violet'
            plt.bar(range(len(regions)), region_site_counts, width=0.4, color=color)
            plt.xticks(range(len(regions)), regions)
            plt.title(f'${dict_mod_display[mod]}$', fontsize=12)
        plt.xlabel('Transcript region', fontsize=12)
        if subplot_ind == 0:
            plt.ylabel('Sites with up- / down-regulation', fontsize=12)
    ymax = np.max(all_region_site_counts)
    ytop = int((ymax // 25 + 1) * 25)
    for subplot_ind in range(len(['psi', 'm6A'])):
        plt.subplot(1, 2, subplot_ind + 1)
        plt.ylim([-ytop, ytop])

    return fig