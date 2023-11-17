import os
import sys
sys.path.append('~/git/mAFiA_dev/mAFiA')
from glob import glob
import pickle
from mAFiA.feature_classifiers import MotifClassifier
from mAFiA.data_containers import Nucleotide

clf_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/ISA-WUE'
out_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/ISA-WUE_new'
os.makedirs(out_dir, exist_ok=True)

for clf_path in glob(os.path.join(clf_dir, '*.pkl')):
    with open(clf_path, 'rb') as h_clf:
        old_clf = pickle.load(h_clf)
    new_clf = MotifClassifier(old_clf.motif, old_clf.classifier_type, old_clf.scaler)
    new_clf.binary_model = old_clf.binary_model
    new_clf.train_nts = [Nucleotide(**old_train_nt.__dict__) for old_train_nt in old_clf.train_nts]
    new_clf.test_nts = [Nucleotide(**old_test_nt.__dict__) for old_test_nt in old_clf.test_nts]
    for k in ['precision', 'recall', 'thresholds', 'auc']:
        new_clf.__dict__[k] = old_clf.__dict__[k]

    new_clf.save(os.path.join(out_dir, f'{new_clf.motif}.pkl'))