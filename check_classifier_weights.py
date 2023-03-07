import os
from joblib import load
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

classifier = 'logistic_regression'
classifier_model_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/MAFIA_classifiers/A_m6A_multiple_noNorm_allReads'
motifs = [
    'GGACA',
    'GGACC',
    'AGACT'
]

classifier_models = {this_motif : load(os.path.join(classifier_model_dir, '{}_{}.joblib'.format(classifier, this_motif))) for this_motif in motifs}

plt.figure(figsize=(16, 9))
for ind, (this_motif, this_model) in enumerate(classifier_models.items()):
    plt.subplot(3, 1, ind+1)
    plt.plot(this_model.coef_.flatten(), label=this_motif)
    plt.xlim([0, len(this_model.coef_.flatten())])
    if ind<(len(motifs)-1):
        plt.xticks([])
    else:
        plt.xlabel('Model Weight', fontsize=15)
    plt.legend(loc='upper left', fontsize=10)