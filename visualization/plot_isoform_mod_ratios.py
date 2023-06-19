import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dataset = 'P2_WT'
annotated_res = f'/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_train_ISA-WUE_test_{dataset}.tsv.merged.annotated'
img_out = f'/home/adrian/img_out/MAFIA/{dataset}_isoforms'
os.makedirs(img_out, exist_ok=True)

df_annotated = pd.read_csv(annotated_res, sep='\t')

num_annotated = (df_annotated['annotation']!='unknown').sum()
percent_annotated = num_annotated / len(df_annotated) * 100
print(f'{num_annotated}/{len(df_annotated)} ({percent_annotated:.1f}%) annotated')

min_coverage = 10

isoform_mod_ratios = {}
for this_index in tqdm(df_annotated['index'].unique()):
    sub_df = df_annotated[df_annotated['index']==this_index]
    sub_df = sub_df[sub_df['annotation']!='unknown']

    isoform_read_nums = [(sub_df['annotation']==this_anno).sum() for this_anno in sub_df['annotation'].unique()]
    if np.sum(np.array(isoform_read_nums)>=min_coverage)<2:
        continue

    isoform_mod_ratios[this_index] = []

    for this_anno in sub_df['annotation'].unique():
        df_anno = sub_df[sub_df['annotation']==this_anno]
        if len(df_anno)<min_coverage:
            continue
        anno_mod_ratio = np.mean(df_anno['mod_prob']>=0.5)
        isoform_mod_ratios[this_index].append((this_anno, anno_mod_ratio, len(df_anno)))

    anno_ratios = [tup[1] for tup in isoform_mod_ratios[this_index]]
    if (np.max(anno_ratios)-np.min(anno_ratios))>=0.5:
        plt.figure(figsize=(5, 5))
        for (this_anno, this_mod_ratio, this_num_read) in isoform_mod_ratios[this_index]:
            plt.hist(sub_df[sub_df['annotation']==this_anno]['mod_prob'].values, range=[0, 1], bins=50, label=f"{this_anno}, {this_mod_ratio:.2f}, {this_num_read}")
        plt.legend(loc='upper right')
        plt.title(f"GLORI {sub_df['Ratio'].values[0]:.2f}")
        plt.savefig(os.path.join(img_out, f"site{sub_df['index'].values[0]}.png"), bbox_inches='tight')
        plt.close('all')