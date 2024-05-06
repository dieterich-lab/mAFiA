import os
from glob import glob
import pickle
import numpy as np

# old_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/DRACH5_retrain'
old_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/psU/PSICO_classifiers/GTTCT_retrain'
new_dir = old_dir + '_stripped'
os.makedirs(new_dir, exist_ok=True)

clf_paths = glob(os.path.join(old_dir, '*.pkl'))
for this_clf_path in clf_paths:
    with open(this_clf_path, 'rb') as h:
        this_clf = pickle.load(h)

        this_clf.train_nts = []
        this_clf.test_nts = []
        this_clf.precision = []
        this_clf.recall = []
        this_clf.thresholds = []
        this_clf.auc = -1
        this_clf.save(os.path.join(new_dir, '{}.pkl'.format(this_clf.motif)))