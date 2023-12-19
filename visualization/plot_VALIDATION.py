import os
from glob import glob
import pickle
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import seaborn as sn

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 7
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################


def get_histogram(mod_probs):
    num_bins = 100
    # num_bins = 20
    bin_edges = np.linspace(0, 1, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_counts, _ = np.histogram(mod_probs, bins=bin_edges)
    norm_counts = bin_counts / np.sum(bin_counts)
    return bin_counts, norm_counts, bin_centers


def calc_prc(modProbs):
    thresholds = np.arange(0, 1.0, 0.01)
    prec = []
    rec = []
    for thresh in thresholds:
        tp = (modProbs['MOD'] > thresh).sum()
        fa = (modProbs['UNM'] > thresh).sum()
        all_pos = len(modProbs['MOD'])

        if (tp==0) and (fa==0):
            break

        prec.append(tp / (tp + fa))
        rec.append(tp / all_pos)
    prec.append(1.0)
    rec.append(0.0)
    area = auc(rec[::-1], prec[::-1])
    return prec, rec, area


def get_sample_sizes(nts, in_labels):
    nt_read_ids = [nt.read_id for nt in nts]
    nt_labels = in_labels.loc[nt_read_ids]['label'].values
    num_unm = (nt_labels == 0).sum()
    num_mod = (nt_labels == 1).sum()
    return num_unm, num_mod

# source_data_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/ISA-WUE'
# read_ids_A_path = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/ISA-WUE_A/read_ids_ISA-WUE_A.txt'
# read_ids_m6A_path = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/ISA-WUE_m6A/read_ids_ISA-WUE_m6A.txt'

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/VALIDATION'

train = 'ISA'
validate = 'WUE'

# train = 'ISA_retrain_GAACT_TGACT'
# validate = 'WUE'

read_ids_A_path = os.path.join(source_data_dir, 'read_ids_ISA-WUE_A.txt')
read_ids_m6A_path = os.path.join(source_data_dir, 'read_ids_ISA-WUE_m6A.txt')

img_out = f'/home/adrian/NCOMMS_revision/images/VALIDATION/train_{train}_validate_{validate}'
os.makedirs(img_out, exist_ok=True)

with open(read_ids_A_path, 'r') as h:
    read_ids_A = [ll.strip('\n') for ll in h.readlines()]
with open(read_ids_m6A_path, 'r') as h:
    read_ids_m6A = [ll.strip('\n') for ll in h.readlines()]

labels = [0] * len(read_ids_A) + [1] * len(read_ids_m6A)
df_labels = pd.DataFrame({'read_id' : read_ids_A + read_ids_m6A, 'label' : labels})
df_labels = df_labels.set_index('read_id')

train_clfs = {}
for clf_path in glob(os.path.join(source_data_dir, train, '*.pkl')):
    with open(clf_path, 'rb') as h_clf:
        clf = pickle.load(h_clf)
        train_clfs[clf.motif] = clf

validate_clfs = {}
for clf_path in glob(os.path.join(source_data_dir, validate, '*.pkl')):
    with open(clf_path, 'rb') as h_clf:
        clf = pickle.load(h_clf)
        validate_clfs[clf.motif] = clf

motif_mod_probs = {}
motif_sample_sizes = {}
for motif, train_clf in train_clfs.items():
    validate_clf = validate_clfs[motif]
    motif_mod_probs[motif] = {}
    motif_sample_sizes[motif] = {}

    validate_features = [test_nt.feature for test_nt in validate_clf.test_nts]
    validate_nt_read_ids = [test_nt.read_id for test_nt in validate_clf.test_nts]
    validate_labels = df_labels.loc[validate_nt_read_ids]['label'].values
    validate_mod_probs = train_clf.binary_model.predict_proba(validate_features)[:, 1]
    motif_mod_probs[motif]['UNM'] = validate_mod_probs[validate_labels==0]
    motif_mod_probs[motif]['MOD'] = validate_mod_probs[validate_labels==1]

    motif_sample_sizes[validate_clf.motif]['train'] = {}
    num_train_unm, num_train_mod = get_sample_sizes(validate_clf.train_nts, df_labels)
    motif_sample_sizes[validate_clf.motif]['train']['UNM'] = num_train_unm
    motif_sample_sizes[validate_clf.motif]['train']['MOD'] = num_train_mod

    motif_sample_sizes[validate_clf.motif]['validate'] = {}
    num_validate_unm, num_validate_mod = get_sample_sizes(validate_clf.test_nts, df_labels)
    motif_sample_sizes[validate_clf.motif]['validate']['UNM'] = num_validate_unm
    motif_sample_sizes[validate_clf.motif]['validate']['MOD'] = num_validate_mod

df_sample_sizes = pd.DataFrame([[motif] + [vv for k, v in counts.items() for kk, vv in v.items()] for motif, counts in motif_sample_sizes.items()])
df_sample_sizes.rename(columns=dict(zip(range(5), ['motif', 'train_unm', 'train_mod', 'validate_unm', 'validate_mod'])), inplace=True)
df_sample_sizes.to_csv(os.path.join(img_out, f'{validate}_train_validate_sample_sizes.tsv'), sep='\t', index=False)

########################################################################################################################
### histogram of prob. distribution ####################################################################################
########################################################################################################################
num_rows = 2
num_cols = 3
fig_width = 8*cm
fig_height = fig_width / gr
x_max = 1.0
y_max = 0
xticks = np.round(np.linspace(0, x_max, 3), 3)

label_colors = {
    'UNM': 'b',
    'MOD': 'r'
}

ordered_motifs = ['GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT']

fig_hist = plt.figure(figsize=(fig_width, fig_height))
axes_hist = fig_hist.subplots(num_rows, num_cols).flatten()
for subplot_ind, this_motif in enumerate(ordered_motifs):
    for label in ['UNM', 'MOD']:
        this_real_count, this_norm_count, this_bin_centers = get_histogram(motif_mod_probs[this_motif][label])
        this_num_samples = len(motif_mod_probs[this_motif][label])
        axes_hist[subplot_ind].step(this_bin_centers, this_real_count, color=label_colors[label], label=f'{label} ({this_num_samples})')

    if subplot_ind>=(num_rows-1)*num_cols:
        axes_hist[subplot_ind].set_xticks(xticks)
    else:
        axes_hist[subplot_ind].set_xticks([])

    axes_hist[subplot_ind].set_title('{}'.format(this_motif.replace('T', 'U')), pad=-10)
    axes_hist[subplot_ind].set_xlim([-0.01, 1.01])
    axes_hist[subplot_ind].legend(loc='upper center', handlelength=1)
fig_hist.tight_layout(pad=0.5)
fig_hist.savefig(os.path.join(img_out, f'hist_oligo_modProbs.{FMT}'), **fig_kwargs)


########################################################################################################################
### confusion matrix ###################################################################################################
########################################################################################################################
THRESH_PROB = 0.5
def calculate_confusion_matrix(mod_probs, thresh_prob=THRESH_PROB):
    mat00 = np.sum(mod_probs['UNM'] < thresh_prob)
    mat01 = np.sum(mod_probs['UNM'] >= thresh_prob)
    mat10 = np.sum(mod_probs['MOD'] < thresh_prob)
    mat11 = np.sum(mod_probs['MOD'] >= thresh_prob)
    return pd.DataFrame([[mat00, mat01], [mat10, mat11]], index=['UNM', 'MOD'], columns=['UNM', 'MOD'])

margin = 0.25
xx = np.array([0, 0, 1, 1]) + margin
yy = np.array([0, 1, 0, 1]) + margin
zz = np.zeros(4)
dx = np.ones(4) - margin
dy = np.ones(4) - margin

fig_confusion = plt.figure(figsize=(18*cm, 12*cm))
# axes_confusion = fig_confusion.subplots(num_rows, num_cols, projection='3d').flatten()
for subplot_ind, this_motif in enumerate(ordered_motifs):
    this_motif_confusion_mat = calculate_confusion_matrix(motif_mod_probs[this_motif])
    dz = this_motif_confusion_mat.values.flatten()
    this_ax = fig_confusion.add_subplot(num_rows, num_cols, subplot_ind+1, projection='3d')
    this_ax.bar3d(xx, yy, zz, dx, dy, dz)
    this_ax.set_xticks([0.5, 1.5], ['UNM', 'MOD'])
    this_ax.xaxis.set_rotate_label(False)
    this_ax.set_xlabel('GT', labelpad=0)
    this_ax.set_yticks([0.5, 1.5], ['UNM', 'MOD'])
    this_ax.yaxis.set_rotate_label(False)
    this_ax.set_ylabel('PRED', labelpad=0)

    for (this_x, this_y, this_z) in zip(xx+0.3, yy+0.15, dz):
        this_ax.text(this_x, this_y, this_z, this_z, c='w')

    # axes_confusion[subplot_ind].matshow(this_motif_confusion_mat)
    # sn.heatmap(ax=axes_confusion[subplot_ind], data=this_motif_confusion_mat, annot=True, fmt='d', cbar=False)
    this_ax.set_title('{}'.format(this_motif.replace('T', 'U')))

fig_confusion.savefig(os.path.join(img_out, f'confusion_matrix_thresh{THRESH_PROB}.{FMT}'), **fig_kwargs)
plt.close('all')
########################################################################################################################
### precision-recall curve #############################################################################################
########################################################################################################################
motif_auprcs = {}

fig_prc = plt.figure(figsize=(fig_width, fig_height))
axes_prc = fig_prc.subplots(num_rows, num_cols).flatten()
for subplot_ind, this_motif in enumerate(ordered_motifs):
    precisions, recalls, auprc = calc_prc(motif_mod_probs[this_motif])
    axes_prc[subplot_ind].plot(recalls, precisions, label='AUC {:.2f}'.format(auprc))
    axes_prc[subplot_ind].set_title('{}'.format(this_motif.replace('T', 'U')), pad=-10)
    axes_prc[subplot_ind].set_xlim([-0.01, 1.01])
    axes_prc[subplot_ind].set_ylim([-0.01, 1.01])
    axes_prc[subplot_ind].legend(loc='lower left', handlelength=1)

    motif_auprcs[this_motif] = auprc

    x_max = 1.0
    y_max = 1.0
    xticks = np.round(np.linspace(0, x_max, 3), 3)
    yticks = np.round(np.linspace(0, y_max, 3), 3)
    if subplot_ind>=(num_rows-1)*num_cols:
        axes_prc[subplot_ind].set_xticks(xticks)
    else:
        axes_prc[subplot_ind].set_xticks([])
    axes_prc[subplot_ind].set_yticks(yticks)
fig_prc.tight_layout(pad=0.5)
fig_prc.savefig(os.path.join(img_out, f'prc_oligo_modProbs.{FMT}'), **fig_kwargs)

motif_auprcs['avg'] = np.mean(list(motif_auprcs.values()))
with open(os.path.join(img_out, 'motif_aucs.tsv'), 'w') as h:
    for this_motif in ordered_motifs:
        h.write(f'{this_motif}\t{np.round(motif_auprcs[this_motif], 3)}\n')
    h.write(f"avg\t{np.round(motif_auprcs['avg'], 3)}\n")