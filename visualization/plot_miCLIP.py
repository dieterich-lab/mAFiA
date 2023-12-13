import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from tqdm import tqdm
from collections import Counter
import pysam
from sklearn.metrics import auc

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/MICLIP'
img_out = '/home/adrian/NCOMMS_revision/images/MICLIP'
os.makedirs(img_out, exist_ok=True)

bam_file = os.path.join(source_data_dir, 'merged.mAFiA.reads.bam')

# miclip_path = '/home/adrian/Data/DRACH/miCLIP_union_flat_exclude_Y_chromosome.ref5mers.bed'
miclip_path = os.path.join(source_data_dir, 'miclip_sites_cov10_6motifs.tsv')

mAFiA_path = os.path.join(source_data_dir, 'merged.mAFiA.sites.bed')
m6Anet_path = os.path.join(source_data_dir, 'data.site_proba.csv.genome.6motifs.combined')
CHEUI_path = os.path.join(source_data_dir, 'site_level_m6A_predictions.txt.genome.6motifs.combined')

# def import_mAFiA(thresh_coverage=10):
#     # dfs = []
#     # for this_chr in chrs:
#     #     dfs.append(pd.read_csv(os.path.join(source_data_dir, f'chr{this_chr}', 'mAFiA.sites.bed'), dtype={'chrom': str}, sep='\t'))
#     # df_mAFiA = pd.concat(dfs, ignore_index=True)
#     df_mAFiA = pd.read_csv(, dtype={'chrom': str}, sep='\t')
#     df_mAFiA_thresh = df_mAFiA[df_mAFiA['coverage']>=thresh_coverage]
#     return df_mAFiA_thresh
#
# def import_m6Anet(thresh_coverage=10):
#     df_m6Anet = pd.read_csv(, dtype={'chrom': str}, sep='\t')
#     df_m6Anet_thresh = df_m6Anet[
#         (df_m6Anet['coverage']>=thresh_coverage)
#         * (df_m6Anet['chrom'].isin([str(x) for x in range(1, 23)] + ['X']))
#         ]
#
#     return df_m6Anet_thresh


def import_results(filepath, thresh_coverage=10):
    df = pd.read_csv(filepath, dtype={'chrom': str}, sep='\t')
    df_thresh = df[
        (df['coverage']>=thresh_coverage)
        * (df['chrom'].isin([str(x) for x in range(1, 23)] + ['X']))
        ]
    return df_thresh

def import_miclip(thresh_source_counts=1):
    df_miclip = pd.read_csv(miclip_path, sep='\t',
                            # names=['chrom', 'chromStart', 'chromEnd', 'source', 'score', 'strand'],
                            dtype={'chrom': str})
    df_miclip = df_miclip[df_miclip['chrom'].isin([str(x) for x in range(1, 23)] + ['X'])]
    df_miclip['modRatio'] = 100
    # df_miclip['source_counts'] = [len(this_source.split(',')) for this_source in df_miclip['source'].values]
    # df_miclip_thresh = df_miclip[df_miclip['source_counts']>=thresh_source_counts]
    return df_miclip


def select_miclip_sites_with_coverage(df_in, thresh_coverage=10):
    with pysam.AlignmentFile(bam_file, "rb" ) as bam:
        coverage = []
        for _, row in tqdm(df_in.iterrows()):
            this_row_coverage = 0
            for pileupcolumn in bam.pileup(row['chrom'], row['chromStart'], row['chromStart']+1, truncate=True):
                if pileupcolumn.pos==row['chromStart']:
                    this_row_coverage = pileupcolumn.n
                    break
            coverage.append(this_row_coverage)
    coverage = np.array(coverage)
    df_in['coverage'] = coverage
    return df_in[df_in['coverage']>=thresh_coverage]

def draw_venn_diagram(dict_dfs, thresh_stoichiometry=50):
    name_sites = {}
    for name, df in dict_dfs.items():
        # df_sel = df[(df['modRatio']>=thresh_stoichiometry)*(df['chrom'].isin(chrs))]
        df_sel = df[(df['modRatio']>=thresh_stoichiometry)]
        name_sites[name] = set([tuple(val) for val in df_sel[['chrom', 'chromStart']].values])

        # num_intersection = len(set(df1_sites).intersection(set(df2_sites)))
        # print(len(df1_sites), len(df2_sites), num_intersection)

    if len(dict_dfs)>=3:
        return venn3(name_sites.values(), name_sites.keys())
    else:
        return venn2(name_sites.values(), name_sites.keys())

df_mAFiA = import_results(mAFiA_path)
df_m6Anet = import_results(m6Anet_path)
df_CHEUI = import_results(CHEUI_path)
df_miclip = import_miclip()

########################################################################################################################
### venn diagram #######################################################################################################
########################################################################################################################
thresh_coverage = 50
sel_motifs = [
    'GGACT',
    'GAACT',
    'GGACA',
    'AGACT',
    'GGACC',
    'TGACT'
]

df_mAFiA_thresh = df_mAFiA[
    df_mAFiA['ref5mer'].isin(sel_motifs)
    * df_mAFiA['coverage']>=thresh_coverage
]

df_m6Anet_thresh = df_m6Anet[
    df_m6Anet['ref5mer'].isin(sel_motifs)
    * df_m6Anet['coverage']>=thresh_coverage
]

df_CHEUI_thresh = df_CHEUI[
    df_CHEUI['ref5mer'].isin(sel_motifs)
    * df_CHEUI['coverage']>=thresh_coverage
]

df_miclip_thresh = df_miclip[
    df_miclip['ref5mer'].isin(sel_motifs)
    * df_miclip['coverage']>=thresh_coverage
]


venn_dict = {
    'mAFiA': df_mAFiA_thresh,
    'm6Anet': df_m6Anet_thresh,
    # 'CHEUI': df_CHEUI_thresh,
    'miClip': df_miclip_thresh,
}

thresh_stoichio = 10
fig_venn, ax = plt.subplots(nrows=1, ncols=1, figsize=(4*cm, 4*cm))
v = draw_venn_diagram(venn_dict, thresh_stoichiometry=thresh_stoichio)
fig_venn.savefig(os.path.join(img_out, f'venn_diagram_cov{thresh_coverage}_stoichio{thresh_stoichio}.{FMT}'), **fig_kwargs)

########################################################################################################################
### precision-recall curve #############################################################################################
########################################################################################################################
def calc_precision_recall_curve(df_gt, df_pred):
    thresh_S = np.linspace(0, 100, 101)
    total_gt = len(df_gt)
    precision = []
    recall = []
    for this_thresh in thresh_S:
        thresh_df_pred = df_pred[df_pred['modRatio']>=this_thresh]
        df_merged = pd.merge(df_gt, thresh_df_pred, on=['chrom', 'chromStart', 'chromEnd'])
        tp = len(df_merged)
        fa = len(thresh_df_pred) - tp
        if (tp==0) and (fa==0):
            precision.append(1.0)
            recall.append(0.0)
        else:
            precision.append(tp / (tp + fa))
            recall.append(tp / total_gt)
    auprc = auc(recall[::-1], precision[::-1])
    return thresh_S, precision, recall, auprc


thresh_mAFiA, prec_mAFiA, recall_mAFiA, auc_mAFiA = calc_precision_recall_curve(df_miclip_thresh, df_mAFiA_thresh)
thresh_m6Anet, prec_m6Anet, recall_m6Anet, auc_m6Anet = calc_precision_recall_curve(df_miclip_thresh, df_m6Anet_thresh)
thresh_CHEUI, prec_CHEUI, recall_CHEUI, auc_CHEUI = calc_precision_recall_curve(df_miclip_thresh, df_CHEUI_thresh)

fig_miclip = plt.figure(figsize=(5*cm, 4*cm))
plt.plot(recall_mAFiA, prec_mAFiA, label=f'mAFiA ({auc_mAFiA:.2f})')
plt.plot(recall_m6Anet, prec_m6Anet, label=f'm6Anet ({auc_m6Anet:.2f})')
plt.plot(recall_CHEUI, prec_CHEUI, label=f'CHEUI ({auc_CHEUI:.2f})')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.xticks(np.linspace(0, 1, 5))
plt.yticks(np.linspace(0, 1, 5))
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend(loc='upper right', handlelength=1)
plt.title(f'Overlap with {len(df_miclip_thresh)} miCLIP sites (cov$\geq${thresh_coverage})')
fig_miclip.savefig(os.path.join(img_out, f'prc_miclip_cov{thresh_coverage}.{FMT}'), **fig_kwargs)
