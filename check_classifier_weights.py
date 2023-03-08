import os
from joblib import load
import pickle
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

NUM_TOP_FEATURES = 10

classifier = 'logistic_regression'
classifier_model_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/MAFIA_classifiers/A_m6A_multiple_noNorm_allReads'
motifs = [
    'GGACA',
    'GGACC',
    'AGACT'
]

classifier_models = {this_motif : load(os.path.join(classifier_model_dir, '{}_{}.joblib'.format(classifier, this_motif))) for this_motif in motifs}
top_model_weights = {}

plt.figure(figsize=(16, 9))
for ind, (this_motif, this_model) in enumerate(classifier_models.items()):
    this_model_weights = this_model.coef_.flatten()
    top_feature_indices = np.argpartition(np.abs(this_model_weights), -NUM_TOP_FEATURES)[-NUM_TOP_FEATURES:]
    top_feature_indices.sort()
    top_model_weights[this_motif] = top_feature_indices
    plt.subplot(3, 1, ind+1)
    plt.plot(this_model_weights, label=this_motif)
    plt.plot(top_feature_indices, this_model_weights[top_feature_indices], 'ro')
    plt.xlim([0, len(this_model_weights)])
    if ind<(len(motifs)-1):
        plt.xticks([])
    else:
        plt.xlabel('Model Weight', fontsize=15)
    # plt.legend(loc='upper left', fontsize=10)
    plt.title('{}: '.format(this_motif) + ', '.join([str(feat_ind) for feat_ind in top_feature_indices]), fontsize=15)
    plt.suptitle('Top {} model weights'.format(NUM_TOP_FEATURES), fontsize=20)
plt.savefig(os.path.join(classifier_model_dir, 'top_{}_model_weights.png'.format(NUM_TOP_FEATURES)), bbox_inches='tight')
plt.close()

with open(os.path.join(classifier_model_dir, 'top_{}_model_weights.pkl'.format(NUM_TOP_FEATURES)), 'wb') as p_out:
    pickle.dump(top_model_weights, p_out)