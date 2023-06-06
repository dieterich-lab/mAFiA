import os
HOME = os.path.expanduser('~')
from glob import glob
import numpy as np
import random
random.seed(0)
from random import sample
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, auc
# import matplotlib
# matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import pickle

class Motif_Classifier:
    def __init__(self, motif, classifier_type, scaler):
        self.motif = motif
        self.classifier_type = classifier_type
        self.scaler = scaler
        if classifier_type == 'svm':
            clf = SVC(gamma='auto', random_state=0, max_iter=1000)
        elif classifier_type == 'logistic_regression':
            clf = LogisticRegression(random_state=0, max_iter=1000)
        else:
            clf = None

        if scaler == 'MaxAbs':
            self.binary_model = make_pipeline(MaxAbsScaler(), clf)
        elif scaler == 'Standard':
            self.binary_model = make_pipeline(StandardScaler(), clf)
        else:
            self.binary_model = clf

    def train(self, unm_nts, mod_nts, frac_test_split=0.25):
        max_num_features = min(len(unm_nts), len(mod_nts))
        sample_unm_nts = sample(unm_nts, max_num_features)
        sample_mod_nts = sample(mod_nts, max_num_features)

        labels = np.array([0 for ii in range(len(sample_unm_nts))] + [1 for ii in range(len(sample_mod_nts))])
        unm_features = [nt.feature for nt in sample_unm_nts]
        mod_features = [nt.feature for nt in sample_mod_nts]
        all_features = np.array(unm_features + mod_features)

        X_train, X_test, y_train, y_test, train_nts, test_nts = train_test_split(all_features, labels, sample_unm_nts+sample_mod_nts, test_size=frac_test_split)
        self.train_nts = train_nts
        self.test_nts = test_nts

        self.binary_model = self.binary_model.fit(X_train, y_train)
        y_score = self.binary_model.decision_function(X_test)
        self.precision, self.recall, self.thresholds = precision_recall_curve(y_test, y_score)
        self.auc = auc(self.recall, self.precision)

        print('AUC {:.2f}'.format(self.auc))

    def test(self, test_nts, mod_thresh=0.5):
        print('Testing {} NTs...'.format(len(test_nts)))
        test_features = [nt.feature for nt in test_nts]
        mod_probs = self.binary_model.predict_proba(test_features)[:, 1]
        for this_nt, this_mod_prob in zip(test_nts, mod_probs):
            this_nt.mod_prob = this_mod_prob
        predictions = np.int32(mod_probs > mod_thresh)
        avg_mod_ratio = np.mean(predictions)
        print('Predicted mod. ratio {:.2f}'.format(avg_mod_ratio))
        return avg_mod_ratio

    def save(self, out_model_path, draw_prc=False):
        with open(out_model_path, 'wb') as h_out:
            pickle.dump(self, h_out, pickle.HIGHEST_PROTOCOL)

        if draw_prc:
            out_img_path = out_model_path.replace('.pkl', '.png')
            plt.figure(figsize=(5, 5))
            plt.plot(self.recall, self.precision, '-')
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.ylim([0, 1.05])
            plt.title('{}\n{} train NTs, {} test NTs\n AUC = {:.2f}'.format(
                self.motif,
                len(self.train_nts),
                len(self.test_nts),
                self.auc
            ))
            plt.savefig(out_img_path, bbox_inches='tight')
            plt.close('all')


def load_motif_classifiers(classifier_dir):
    print('Loading motif classifiers...')
    classifier_paths = glob(os.path.join(classifier_dir, '*.pkl'))
    motif_classifiers = {}
    for this_path in classifier_paths:
        with open(this_path, 'rb') as h_in:
            this_motif_classifier = pickle.load(h_in)
        motif_classifiers[this_motif_classifier.motif] = this_motif_classifier
    classifier_motifs = list(motif_classifiers.keys())
    print('Target motifs: {}'.format(', '.join(classifier_motifs)))

    return motif_classifiers