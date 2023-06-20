import os
import pandas as pd
import numpy as np
from tqdm import tqdm

#######################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt

cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
# FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

# dataset = 'P2_WT'
dataset = '100_WT_0_IVT'
annotated_res = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{dataset}.tsv.merged.annotated'
img_out = f'/home/adrian/img_out/MAFIA/{dataset}_isoforms_{FMT}'
os.makedirs(img_out, exist_ok=True)

df_annotated = pd.read_csv(annotated_res, sep='\t')

num_annotated = (df_annotated['annotation']!='unknown').sum()
percent_annotated = num_annotated / len(df_annotated) * 100
print(f'{num_annotated}/{len(df_annotated)} ({percent_annotated:.1f}%) annotated')

thresh_p_val = 1E-99
min_coverage_site = 50
min_coverage_isoform = 10
gap_ratio = 0.1

isoform_mod_ratios = {}
for this_index in tqdm(df_annotated['index'].unique()):
    sub_df = df_annotated[df_annotated['index']==this_index]
    sub_df = sub_df[sub_df['annotation']!='unknown']
    if len(sub_df)<min_coverage_site:
        continue

    isoform_read_nums = [(sub_df['annotation']==this_anno).sum() for this_anno in sub_df['annotation'].unique()]
    if (float(sub_df['Pvalue'].values[0])>thresh_p_val) or (np.sum(np.array(isoform_read_nums)>=min_coverage_isoform)<2):
        continue

    isoform_mod_ratios[this_index] = []

    for this_anno in sub_df['annotation'].unique():
        df_anno = sub_df[sub_df['annotation']==this_anno]
        if len(df_anno)<min_coverage_isoform:
            continue
        anno_mod_ratio = np.mean(np.float64(df_anno['mod_prob'].values)>=gap_ratio)
        isoform_mod_ratios[this_index].append((this_anno, anno_mod_ratio, len(df_anno)))

    anno_ratios = [tup[1] for tup in isoform_mod_ratios[this_index]]
    if (np.max(anno_ratios)-np.min(anno_ratios))>=gap_ratio:
        chr = sub_df['Chr'].unique()[0]
        pos = sub_df['Sites'].unique()[0]
        gene = sub_df['Gene'].unique()[0]
        glori_ratio = float(sub_df['Ratio'].values[0])

        ### hist. mod. prob. ###
        # plt.figure(figsize=(5*cm, 5*cm))
        # for (this_anno, this_mod_ratio, this_num_read) in isoform_mod_ratios[this_index]:
        #     plt.hist(np.float64(sub_df[sub_df['annotation']==this_anno]['mod_prob'].values), range=[0, 1], bins=20, alpha=0.5, label=f"{this_anno}")
        # plt.xlabel('Mod. Prob.')
        # plt.ylabel('Read Counts')
        # plt.legend(loc='upper left')
        # plt.title(f"{chr}, pos {pos}, gene {gene}\nGLORI {float(sub_df['Ratio'].values[0]):.2f}, mAFiA {float(sub_df['pred_mod_ratio'].values[0]):.2f}")
        # plt.savefig(os.path.join(img_out, f"hist_site{sub_df['index'].values[0]}.{FMT}"), **fig_kwargs)

        ### dot plot ###
        def set_row_xs(num_dots, x_spacing=0.08):
            if num_dots % 2 == 0:
                half_width = (num_dots / 2 - 0.5) * x_spacing
            else:
                half_width = (num_dots // 2) * x_spacing
            return np.linspace(-half_width, half_width, num_dots)

        xrange = [-1, 1]
        row_max = 50
        default_x_spacing = (xrange[1] - xrange[0]) / row_max
        y_delta = 0.005
        y_round_up = 0.1
        dot_size = 0.5

        cmap = plt.cm.get_cmap('rainbow', len(isoform_mod_ratios[this_index]))

        fig = plt.figure(figsize=(4*cm, 6*cm))
        ax1 = fig.add_subplot()
        for iso_ind, (this_anno, this_mod_ratio, this_num_read) in enumerate(isoform_mod_ratios[this_index]):
            if this_num_read<=row_max:
                xs = set_row_xs(this_num_read, default_x_spacing)
                ys = [this_mod_ratio] * this_num_read
            else:
                xs = []
                ys = []
                row_dots = row_max
                row_y = this_mod_ratio
                while this_num_read>0:
                    xs.append(set_row_xs(row_dots, default_x_spacing))
                    ys.append([row_y] * row_dots)
                    this_num_read -= row_dots
                    row_dots = min(row_dots-1, this_num_read)
                    row_y += y_delta
                xs = np.concatenate(xs)
                ys = np.concatenate(ys)
            ax1.scatter(xs, ys, color=cmap(iso_ind), s=dot_size)

        mafia_weighted = np.sum([tup[1] * tup[2] for tup in isoform_mod_ratios[this_index]]) / np.sum(
            [tup[2] for tup in isoform_mod_ratios[this_index]])

        ytick_pos = [tup[1] for tup in isoform_mod_ratios[this_index]] + [glori_ratio]
        yticks = [tup[0] for tup in isoform_mod_ratios[this_index]] + ['GLORI']
        ylim = [(min(ytick_pos)//y_round_up)*y_round_up - 0.05, (max(ytick_pos)//y_round_up+1)*y_round_up + 0.05]
        ax1.set_xticks([])
        ax1.set_yticks(ytick_pos, yticks)

        ax1.axhline(glori_ratio, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        ax1.set_title(f"{chr}: pos {pos}\ngene {gene}")
        ax1.set_ylim(ylim)
        ax2 = ax1.twinx()
        ax2.set_yticks(ytick_pos, [f"{ytick:.2f}" for ytick in ytick_pos])
        ax2.set_ylim(ylim)
        ax2.set_ylabel('Mod. Ratio', rotation=-90, labelpad=10)

        plt.savefig(os.path.join(img_out, f"dot_site{sub_df['index'].values[0]}.{FMT}"), **fig_kwargs)
        plt.close('all')

### histogram of splitting amount ###
splitting_amount = []
for k, v in isoform_mod_ratios.items():
    if len(v)>1:
        mod_ratios = [tup[1] for tup in v]
        splitting_amount.append(np.max(mod_ratios) - np.min(mod_ratios))

plt.figure(figsize=(5*cm, 5*cm))
plt.hist(splitting_amount, range=[0, 1], bins=50)
plt.xlabel('Max. difference between isoform mod. ratios')
plt.ylabel('Site counts')
plt.xticks([0, 0.5, 1.0])
plt.yticks([0, 100, 200])
plt.savefig(os.path.join(img_out, f"hist_max_diff_isoform_mod_ratios.{FMT}"), **fig_kwargs)
plt.close()

### dominant isoform to GLORI ###
# mod_ratio_pairs = []
# for k, v in isoform_mod_ratios.items():
#     glori_ratio = float(df_annotated[df_annotated['index']==k]['Ratio'].values[0])
#     dominant_mod_ratio = v[np.argmax([tup[2] for tup in v])][1]
#     mod_ratio_pairs.append((glori_ratio, dominant_mod_ratio))
#
# plt.figure(figsize=(5*cm, 5*cm))
# plt.scatter([pair[0] for pair in mod_ratio_pairs], [pair[1] for pair in mod_ratio_pairs])
