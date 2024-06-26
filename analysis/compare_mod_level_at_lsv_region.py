import os
import pandas as pd
pd.set_option('display.max_columns', None)
from pybedtools import BedTool
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

# gtf_file = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf'
# gtf = BedTool(gtf_file)

thresh_conf = 0.0
# day = 'day1'

# majiq_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/majiq/TAC/voila_modulize_{}/summary.tsv'
sham_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_{}/chrALL.mAFiA.sites.bed'
tac_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_{}/chrALL.mAFiA.sites.bed'
img_out = '/home/adrian/img_out/TAC_splice_site_mod_changes'
os.makedirs(img_out, exist_ok=True)

# df_majiq = pd.read_csv(majiq_file, sep='\t', comment='#')
# df_sham = pd.read_csv(sham_file, sep='\t', dtype={'chrom': str})
# df_tac = pd.read_csv(tac_file, sep='\t', dtype={'chrom': str})

def filter_df(in_df, chrom, strand, chromStart, chromEnd, thresh_conf):
    return in_df[
        (in_df['chrom'] == chrom)
        * (in_df['strand'] == strand)
        * (in_df['chromStart'] >= chromStart)
        * (in_df['chromEnd'] <= chromEnd)
        * (in_df['confidence'] >= thresh_conf)
    ]

mods = ['m6A', 'psi']

gene_exon_ranges = {
    'Fhl1': {
        '3': [56754335, 56754391],
        # '4': [56779437, 56779523],
        '5': [56779723, 56779794],
        # '6': [56779819, 56779908],
        # '7': [56786527, 56786582],
        # '8': [56787725, 56787836],
        '9': [56787989, 56788170],
        # 'intron': [56788171, 56789143],
        # '10': [56789144, 56789743]
    },
    'Synpo2l': {
        '3': [20665744, 20666258],
        # '3intron4': [20664575, 20665743],
        '4': [20664397, 20664574],
        '4intron5': [20662466, 20664396],
        '5': [20658946, 20662465]
    },
    'Rcan1': {
        '2': [92465833, 92466146],
        '4': [92399896, 92400077],
        '6': [92397263, 92397436]
    }
}

# for _, row in df_majiq.iterrows():
#     gene_name = row['gene_name']
    # if gene_name not in ['Fhl1', 'Synpo2l']:
    # if gene_name not in ['Fhl1']:
    #     continue

    # exon_ranges = gene_exon_ranges[gene_name]
    # chrom = row['seqid']
    # strand = row['strand']
    # lsv_regions = row['lsv_id'].split(';')
    # for this_lsv_region in lsv_regions:
    #     source_target = this_lsv_region.split(':')[-2]
    #     chromStart, chromEnd = [int(x) for x in this_lsv_region.split(':')[-1].split('-')]

# gene_name = 'Fhl1'
# exon_ranges = gene_exon_ranges[gene_name]
# chrom = 'X'
# strand = '+'

# gene_name = 'Synpo2l'
# exon_ranges = gene_exon_ranges[gene_name]
# chrom = '14'
# strand = '-'

gene_name = 'Rcan1'
exon_ranges = gene_exon_ranges[gene_name]
chrom = '16'
strand = '-'

sel_days = ['day1', 'day7', 'day21', 'day56']
for this_exon_ind, this_exon_range in exon_ranges.items():
    chromStart, chromEnd = this_exon_range
    chromStart -= 1
    chromEnd -= 1
    df_days = {}
    for day in sel_days:
        df_sham = pd.read_csv(sham_file.format(day), sep='\t', dtype={'chrom': str})
        df_tac = pd.read_csv(tac_file.format(day), sep='\t', dtype={'chrom': str})
        sub_df_sham = filter_df(df_sham, chrom, strand, chromStart, chromEnd, thresh_conf)
        sub_df_tac = filter_df(df_tac, chrom, strand, chromStart, chromEnd, thresh_conf)
        merged_df_sham_tac = pd.merge(sub_df_sham, sub_df_tac,
                                      on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                                      suffixes=['_sham', '_tac'],
                                      how='outer'
                                      )
        merged_df_sham_tac.rename(inplace=True, columns={
            'modRatio_sham': f'modRatio_sham_{day}',
            'confidence_sham': f'confidence_sham_{day}',
            'coverage_sham': f'coverage_sham_{day}',
            'modRatio_tac': f'modRatio_tac_{day}',
            'confidence_tac': f'confidence_tac_{day}',
            'coverage_tac': f'coverage_tac_{day}'
        })
        df_days[day] = merged_df_sham_tac

        num_rows = len(mods)
        num_cols = 1
        merged_df_sham_tac[f'delta_{day}'] = merged_df_sham_tac[f'modRatio_tac_{day}'] - merged_df_sham_tac[f'modRatio_sham_{day}']

        if len(merged_df_sham_tac)>0:
            plt.figure(figsize=(5 * num_cols, 4 * num_rows))
            for this_row_ind, this_mod in enumerate(mods):
                plt.subplot(num_rows, num_cols, this_row_ind + 1)
                this_df = merged_df_sham_tac[merged_df_sham_tac['name'] == this_mod]
                labelled = False
                for _, this_site in this_df.iterrows():
                    x = this_site['chromEnd']
                    if not labelled:
                        plt.plot(x, this_site[f'modRatio_sham_{day}'], 'b+', label=f'SHAM_{day}')
                        plt.plot(x, this_site[f'modRatio_tac_{day}'], 'r+', label=f'TAC_{day}')
                        labelled = True
                    else:
                        plt.plot(x, this_site[f'modRatio_sham_{day}'], 'b+')
                        plt.plot(x, this_site[f'modRatio_tac_{day}'], 'r+')
                    lc = 'r' if this_site[f'delta_{day}'] >= 0 else 'b'
                    plt.plot((x, x), (this_site[f'modRatio_sham_{day}'], this_site[f'modRatio_tac_{day}']), '-', c=lc)

                plt.ylim([-5, 100])
                if this_row_ind == (num_rows - 1):
                    plt.xlabel('Genome coord.', fontsize=12)
                plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
                plt.axvline(chromStart, c='gray')
                plt.axvline(chromEnd, c='gray')
                if strand == '+':
                    plt.xticks([chromStart, chromEnd],
                               [f"5'\nchr{chrom}: {chromStart + 1}", f"3'\nchr{chrom}: {chromEnd + 1}"])
                    legend_loc = 'upper left'
                else:
                    plt.xticks([chromStart, chromEnd],
                               [f"3'\nchr{chrom}: {chromStart + 1}", f"5'\nchr{chrom}: {chromEnd + 1}"])
                    legend_loc = 'upper right'
                if this_row_ind == 0:
                    plt.legend(loc=legend_loc)
                region_len = chromEnd - chromStart
                plt.xlim([chromStart - 0.1 * region_len, chromEnd + 0.1 * region_len])
                if this_row_ind == 0:
                    plt.title(f"exon {this_exon_ind}", fontsize=12)
            plt.suptitle(f'{gene_name}\n{day}', fontsize=15)

            plt.savefig(os.path.join(img_out, f'{gene_name}_{day}.png'), bbox_inches='tight')
            plt.close('all')

    df_all_days = pd.merge(df_days['day1'], df_days['day7'],
                           on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                           )
    df_all_days = pd.merge(df_all_days, df_days['day21'],
                           on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                           )
    df_all_days = pd.merge(df_all_days, df_days['day56'],
                           on=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'ref5mer'],
                           )

    for _, this_row in df_all_days.iterrows():
        plt.figure(figsize=(5, 4))
        plt.plot(range(1, len(sel_days)+1), [this_row['modRatio_sham_day1'], this_row['modRatio_sham_day7'], this_row['modRatio_sham_day21'], this_row['modRatio_sham_day56']], 'b-o', label='SHAM')
        plt.plot(range(1, len(sel_days)+1), [this_row['modRatio_tac_day1'], this_row['modRatio_tac_day7'], this_row['modRatio_tac_day21'], this_row['modRatio_tac_day56']], 'r-o', label='TAC')
        plt.legend(loc='upper left')
        plt.xticks(range(1, len(sel_days)+1), sel_days)
        plt.xlabel('Day', fontsize=12)
        plt.ylabel('S', fontsize=12)
        plt.title(f"chr{this_row['chrom']}: {this_row['chromEnd']}, ${{{dict_mod_display[this_row['name']]}}}$", fontsize=15)
        plt.savefig(
            os.path.join(img_out,
                         f"{gene_name}_chr{this_row['chrom']}_{this_row['chromEnd']}_{this_row['name']}.png"),
            bbox_inches='tight'
        )
        plt.close('all')

    cond_colors = {
        'tac': 'r',
        'sham': 'b'
    }

    labels = [(mpatches.Patch(color=cond_colors[this_label]), this_label.upper()) for this_label in ['tac', 'sham']]

    plt.figure(figsize=(8, 8))
    for subplot_ind, this_mod in enumerate(['m6A', 'psi']):
        plt.subplot(2, 1, subplot_ind+1)
        for cond_ind, condition in enumerate(['tac', 'sham']):
            for day_ind, this_day in enumerate(sel_days):
                label = condition.upper() if day_ind == 0 else None
                this_data = df_all_days[df_all_days['name'] == this_mod][f'modRatio_{condition}_{this_day}'].values
                this_data = this_data[~np.isnan(this_data)]
                if len(this_data):
                    violin_parts = plt.violinplot(this_data, [cond_ind+1-0.5+day_ind*0.25], widths=0.2, showmeans=True)
                    for pc in violin_parts['bodies']:
                        pc.set_facecolor(cond_colors[condition])
                        pc.set_edgecolor(cond_colors[condition])
                    violin_parts['cmaxes'].set_edgecolor(cond_colors[condition])
                    violin_parts['cmeans'].set_edgecolor(cond_colors[condition])
                    violin_parts['cmins'].set_edgecolor(cond_colors[condition])
                    violin_parts['cbars'].set_edgecolor(cond_colors[condition])
        plt.xticks([this_cond_ind+1-0.5+this_day_ind*0.25 for this_cond_ind in range(2) for this_day_ind in range(len(sel_days))], sel_days + sel_days)
        # plt.ylim([-5, 100])
        plt.ylabel(f'$S_{{{dict_mod_display[this_mod]}}}$', fontsize=12)
        plt.xlabel('Days', fontsize=12)
        plt.legend(*zip(*labels))
    plt.suptitle(f'{gene_name}, exon {this_exon_ind}, chr{chrom}: {chromStart+1}-{chromEnd+1}', fontsize=15)
    plt.savefig(os.path.join(img_out, f'mod_dist_{gene_name}_exon{this_exon_ind}.png'), bbox_inches='tight')