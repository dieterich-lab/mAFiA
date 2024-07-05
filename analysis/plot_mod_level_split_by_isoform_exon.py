import os
import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC'
img_out = '/home/adrian/img_out/TAC_mod_level_by_isoform_exon'
os.makedirs(img_out, exist_ok=True)

### Fhl1 ###############################################################################################################
gene_name = 'Fhl1'
chrom = 'X'
strand = '+'
exon_ranges = {
    '1': [56731787, 56731868],
    '3': [56754335, 56754391],
    'new': [56764288, 56764402],
    '5': [56779723, 56779794],
    '9': [56787989, 56788170],
    '10': [56789144, 56789743],
    '11': [56789877, 56790046],
    '12': [56790464, 56790666],
    '13': [56791844, 56793346]
}
# isoforms = ['ENSMUST00000023854', 'ENSMUST00000114772', 'unclassified.new_exon']
isoforms = ['ENSMUST00000023854', 'ENSMUST00000114772']
events = isoforms
# colors = ['red', 'purple', 'green']
colors = ['red', 'purple']

### Rcan1 ##############################################################################################################
# gene_name = 'Rcan1'
# chrom = '16'
# strand = '-'
# exon_ranges = {
#     '2': [92465833, 92466146],
#     '4': [92399896, 92400077],
#     '6': [92397263, 92397436],
#     '7': [92395866, 92396025],
#     '8': [92391953, 92393628]
# }
#
# events = ['exon2', 'exon4']
# isoforms = ['ENSMUST00000060005', 'ENSMUST00000023672']
# colors = ['blue', 'purple']

########################################################################################################################
mods = ['m6A', 'psi']
conditions = ['TAC', 'SHAM']
days = ['day1', 'day7', 'day21', 'day56']
dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}
dict_event_isoform = {ev: iso for (ev, iso) in zip(events, isoforms)}
dict_event_color = {ev: color for (ev, color) in zip(events, colors)}
dict_condition_color = {'TAC': 'red', 'SHAM': 'blue'}

### split by condition and day ###
isoform_labels = [(mpatches.Patch(color=dict_event_color[this_event]), dict_event_isoform[this_event]) for this_event in events]

for this_cond in conditions:
    for this_day in days:
        plt.figure(figsize=(10, 10))
        for mod_ind, this_mod in enumerate(mods):
            plt.subplot(2, 1, mod_ind+1)
            for event_ind, this_event in enumerate(events):
                bed_file = os.path.join(res_dir, '_'.join([this_cond, this_day]), 'bambu', f'{gene_name}.{this_event}.bed')
                if os.path.exists(bed_file):
                    df_bed = pd.read_csv(bed_file, sep='\t', dtype={'chrom': str})
                    df_mod = df_bed[df_bed['name']==this_mod]
                    exon_mod_ratios = [
                        df_mod[(df_mod['chromStart'] >= (this_exon_range[0]-1))
                               * (df_mod['chromEnd'] <= (this_exon_range[1]-1))
                               ]['modRatio'].values
                        for this_exon_range in exon_ranges.values()
                        ]
                    exon_mod_ratios = [this_arr if len(this_arr) else np.array([0]) for this_arr in exon_mod_ratios]
                else:
                    exon_mod_ratios = [np.array([0]) for this_exon_range in exon_ranges]

                violin_parts = plt.violinplot(exon_mod_ratios,
                                              np.arange(len(exon_ranges))-0.1+0.2*event_ind, widths=0.2)
                for pc in violin_parts['bodies']:
                    pc.set_facecolor(dict_event_color[this_event])
                    pc.set_edgecolor(dict_event_color[this_event])
                violin_parts['cmaxes'].set_edgecolor(dict_event_color[this_event])
                # violin_parts['cmeans'].set_edgecolor(dict_event_color[this_event])
                violin_parts['cmins'].set_edgecolor(dict_event_color[this_event])
                violin_parts['cbars'].set_edgecolor(dict_event_color[this_event])
            # plt.legend(loc='upper left', fontsize=10)
            plt.xticks(range(len(exon_ranges)), list(exon_ranges.keys()))
            plt.xlabel('Exon', fontsize=12)
            plt.ylim([-5, 105])
            plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
            # plt.title(f'${{{}}}$', fontsize=15)
            if mod_ind==0:
                plt.legend(*zip(*isoform_labels), loc='upper left')
        plt.suptitle(f'{this_cond} {this_day}', fontsize=15)
        plt.savefig(os.path.join(img_out, f'{gene_name}_{this_cond}_{this_day}.png'))
        plt.close('all')


### split by day and isoform ###
cond_labels = [(mpatches.Patch(color=dict_condition_color[this_cond]), this_cond) for this_cond in conditions]

for this_day in days:
    for event_ind, this_event in enumerate(events):
        plt.figure(figsize=(10, 10))
        for mod_ind, this_mod in enumerate(mods):
            plt.subplot(2, 1, mod_ind+1)
            for cond_ind, this_cond in enumerate(conditions):
                bed_file = os.path.join(res_dir, '_'.join([this_cond, this_day]), 'bambu', f'{gene_name}.{this_event}.bed')
                if os.path.exists(bed_file):
                    df_bed = pd.read_csv(bed_file, sep='\t', dtype={'chrom': str})
                    df_mod = df_bed[df_bed['name']==this_mod]
                    exon_mod_ratios = [
                        df_mod[(df_mod['chromStart'] >= (this_exon_range[0]-1))
                               * (df_mod['chromEnd'] <= (this_exon_range[1]-1))
                               ]['modRatio'].values
                        for this_exon_range in exon_ranges.values()
                        ]
                    exon_mod_ratios = [this_arr if len(this_arr) else np.array([0]) for this_arr in exon_mod_ratios]
                else:
                    exon_mod_ratios = [np.array([0]) for this_exon_range in exon_ranges]

                violin_parts = plt.violinplot(exon_mod_ratios,
                                              np.arange(len(exon_ranges))-0.1+0.2*cond_ind, widths=0.2)
                for pc in violin_parts['bodies']:
                    pc.set_facecolor(dict_condition_color[this_cond])
                    pc.set_edgecolor(dict_condition_color[this_cond])
                violin_parts['cmaxes'].set_edgecolor(dict_condition_color[this_cond])
                # violin_parts['cmeans'].set_edgecolor(dict_event_color[this_event])
                violin_parts['cmins'].set_edgecolor(dict_condition_color[this_cond])
                violin_parts['cbars'].set_edgecolor(dict_condition_color[this_cond])
            # plt.legend(loc='upper left', fontsize=10)
            plt.xticks(range(len(exon_ranges)), list(exon_ranges.keys()))
            plt.xlabel('Exon', fontsize=12)
            plt.ylim([-5, 105])
            plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
            if mod_ind==0:
                plt.title(dict_event_isoform[this_event], fontsize=15)
            if mod_ind==0:
                plt.legend(*zip(*cond_labels), loc='upper left')
        plt.suptitle(this_day, fontsize=15)
        plt.savefig(os.path.join(img_out, f'{gene_name}_{dict_event_isoform[this_event]}_{this_day}.png'))
        plt.close('all')


### split by exon ###
conditions = ['TAC']
for this_exon_ind, this_exon_range in exon_ranges.items():
    plt.figure(figsize=(6*len(events), 10))
    for mod_ind, this_mod in enumerate(mods):
        for event_ind, this_event in enumerate(events):
            plt.subplot(2, len(events), mod_ind * len(events) + event_ind + 1)
            for cond_ind, this_cond in enumerate(conditions):
                day_mod_ratios = []
                for day_ind, this_day in enumerate(days):
                    bed_file = os.path.join(res_dir, '_'.join([this_cond, this_day]), 'bambu', f'{gene_name}.{this_event}.bed')
                    if os.path.exists(bed_file):
                        df_bed = pd.read_csv(bed_file, sep='\t', dtype={'chrom': str})
                        df_mod = df_bed[df_bed['name']==this_mod]
                        exon_mod_ratios = df_mod[
                            (df_mod['chromStart'] >= (this_exon_range[0]-1))
                            * (df_mod['chromEnd'] <= (this_exon_range[1]-1))
                            ]['modRatio'].values
                        exon_mod_ratios = exon_mod_ratios if len(exon_mod_ratios) else np.array([0])
                    else:
                        exon_mod_ratios = np.array([0])
                    day_mod_ratios.append(exon_mod_ratios)

                violin_parts = plt.violinplot(day_mod_ratios,
                                              cond_ind - 0.3 + np.arange(len(days))*0.2, widths=0.1,
                                              # quantiles=[[0.75]] * len(days)
                                              )
                this_color = dict_event_color[this_event]
                for pc in violin_parts['bodies']:
                    pc.set_facecolor(this_color)
                    pc.set_edgecolor(this_color)
                violin_parts['cmaxes'].set_edgecolor(this_color)
                # violin_parts['cquantiles'].set_edgecolor(this_color)
                violin_parts['cmins'].set_edgecolor(this_color)
                violin_parts['cbars'].set_edgecolor(this_color)
            # if (mod_ind==0) and (event_ind==0):
            #     plt.legend(*zip(*cond_labels), loc='upper left')
            plt.xticks(np.concatenate([cond_ind - 0.3 + np.arange(len(days))*0.2 for cond_ind in range(len(conditions))]),
                       days*len(conditions))
            if mod_ind==(len(mods)-1):
                plt.xlabel('Days', fontsize=12)
            plt.ylim([-5, 105])
            plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
            plt.title(f'{dict_event_isoform[this_event]}', fontsize=15)
            # if (mod_ind==0) and (cond_ind==0):
            #     plt.legend(*zip(*labels), loc='lower left')
            plt.suptitle(f'Exon {this_exon_ind}', fontsize=20)
    plt.savefig(os.path.join(img_out, f'{gene_name}_exon{this_exon_ind}.png'))
    plt.close('all')
