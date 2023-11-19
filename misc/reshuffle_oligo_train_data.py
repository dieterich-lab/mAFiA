import os
from glob import glob
import pickle
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, auc
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MaxAbsScaler

img_out = '/home/adrian/img_out/NCOMMS_rev/classifier_validation'

def get_labels():
    read_ids_A_path = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/ISA-WUE_A/read_ids_ISA-WUE_A.txt'
    read_ids_m6A_path = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/ISA-WUE_m6A/read_ids_ISA-WUE_m6A.txt'
    with open(read_ids_A_path, 'r') as h:
        read_ids_A = [ll.strip('\n') for ll in h.readlines()]
    with open(read_ids_m6A_path, 'r') as h:
        read_ids_m6A = [ll.strip('\n') for ll in h.readlines()]
    labels = [0] * len(read_ids_A) + [1] * len(read_ids_m6A)
    df_labels = pd.DataFrame({'read_id': read_ids_A + read_ids_m6A, 'label': labels})
    return df_labels.set_index('read_id')


def get_sample_sizes(nts, df_labels):
    nt_read_ids = [nt.read_id for nt in nts]
    nt_labels = df_labels.loc[nt_read_ids]['label'].values
    num_unm = (nt_labels == 0).sum()
    num_mod = (nt_labels == 1).sum()
    return (num_unm, num_mod)


clf_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/ISA-WUE'
labels = get_labels()

# motif_sample_sizes = {}
# for clf_path in glob(os.path.join(clf_dir, '*.pkl')):
#     with open(clf_path, 'rb') as h_clf:
#         clf = pickle.load(h_clf)
#         print(clf.motif)
#         motif_sample_sizes[clf.motif] = {}
#
#         motif_sample_sizes[clf.motif]['train'] = {}
#         num_train_unm, num_train_mod = get_sample_sizes(clf.train_nts, labels)
#         motif_sample_sizes[clf.motif]['train']['UNM'] = num_train_unm
#         motif_sample_sizes[clf.motif]['train']['MOD'] = num_train_mod
#
#         motif_sample_sizes[clf.motif]['validate'] = {}
#         num_validate_unm, num_validate_mod = get_sample_sizes(clf.test_nts, labels)
#         motif_sample_sizes[clf.motif]['validate']['UNM'] = num_validate_unm
#         motif_sample_sizes[clf.motif]['validate']['MOD'] = num_validate_mod
#
# df_sample_sizes = pd.DataFrame([[motif] + [vv for k, v in counts.items() for kk, vv in v.items()] for motif, counts in motif_sample_sizes.items()])
# df_sample_sizes.rename(columns=dict(zip(range(5), ['motif', 'train_unm', 'train_mod', 'validate_unm', 'validate_mod'])), inplace=True)
# df_sample_sizes.to_csv(os.path.join(img_out, 'oligo_train_validate_sample_sizes.tsv'), sep='\t', index=False)

########################################################################################################################
### reschuffle #########################################################################################################
########################################################################################################################
def get_train_validate_auc(train_nts, validate_nts, df_labels):
    train_features = [train_nt.feature for train_nt in train_nts]
    train_read_ids = [train_nt.read_id for train_nt in train_nts]
    train_labels = df_labels.loc[train_read_ids]['label'].values
    validate_features = [validate_nt.feature for validate_nt in validate_nts]
    validate_read_ids = [validate_nt.read_id for validate_nt in validate_nts]
    validate_labels = df_labels.loc[validate_read_ids]['label'].values

    binary_model = make_pipeline(MaxAbsScaler(), LogisticRegression(random_state=0, max_iter=1000, verbose=True))
    binary_model = binary_model.fit(train_features, train_labels)
    y_score = binary_model.decision_function(validate_features)
    precision, recall, thresholds = precision_recall_curve(validate_labels, y_score)
    return auc(recall, precision)

motif_aucs = {}
for clf_path in glob(os.path.join(clf_dir, '*.pkl')):
    with open(clf_path, 'rb') as h_clf:
        this_clf = pickle.load(h_clf)

    this_motif = this_clf.motif
    motif_aucs[this_motif] = []
    print(this_motif)

    this_trial_train_nts = this_clf.train_nts
    this_trial_validate_nts = this_clf.test_nts
    num_train = len(this_clf.train_nts)
    num_validate = len(this_clf.test_nts)
    this_trial_auc = get_train_validate_auc(this_trial_train_nts, this_trial_validate_nts, labels)
    motif_aucs[this_motif].append(this_trial_auc)
    print(f'Trial 0: AUC {this_trial_auc:.2f}')

    breaks = list(range(0, num_train, num_validate)) + [num_train]
    parts_train_nts = [this_clf.train_nts[breaks[i]:breaks[i+1]] for i in range(3)]

    for this_trial in range(3):
        this_trial_train_nts = this_clf.test_nts + [this_nt for x in range(3) if x!=this_trial for this_nt in parts_train_nts[x]]
        this_trial_validate_nts = parts_train_nts[this_trial]
        this_trial_auc = get_train_validate_auc(this_trial_train_nts, this_trial_validate_nts, labels)
        motif_aucs[this_motif].append(this_trial_auc)
        print(f'Trial {this_trial+1}: AUC {this_trial_auc:.2f}')

df_reshuffle_aucs = pd.DataFrame([[motif] + aucs for motif, aucs in motif_aucs.items()])
df_reshuffle_aucs.rename(columns=dict(zip(range(5), ['motif', 'trial1', 'trial2', 'trial3', 'trial4'])), inplace=True)
df_reshuffle_aucs.to_csv(os.path.join(img_out, 'reshuffle_aucs.tsv'), sep='\t', index=False)