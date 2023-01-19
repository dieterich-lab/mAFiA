import os
HOME = os.path.expanduser('~')
import numpy as np
from scipy.spatial.distance import pdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from collections import Counter
import matplotlib.pyplot as plt

img_out = os.path.join(HOME, 'img_out')

def cluster_features(id_features):
    all_ids = []
    all_features = []
    for k, v in id_features.items():
        all_ids.append(k)
        all_features.append(v)

    ### check features ###
    # plt.figure(figsize=(16, 16))
    # for i in range(16):
    #     plt.subplot(4, 4, i+1)
    #     plt.plot(all_features[i])
    # plt.savefig(os.path.join(img_out, 'feature_vectors.png'), bbox_inches='tight')
    # plt.close('all')

    num_features = len(all_features)
    vec_w = 1.0 - pdist(np.vstack(all_features), metric='cosine')

    # plt.figure(figsize=(8, 8))
    # plt.hist(vec_w, bins=50)
    # plt.savefig(os.path.join(img_out, 'hist_cosine_similarities.png'), bbox_inches='tight')
    # plt.close('all')

    ### cluster by connected components ###
    threshold = 0.5
    vec_i, vec_j = np.triu_indices(num_features, k=1)
    mask = (vec_w>=threshold)
    sel_vec_i = vec_i[mask]
    sel_vec_j = vec_j[mask]
    sel_vec_w = vec_w[mask]

    mat_w = coo_matrix((sel_vec_w, (sel_vec_i, sel_vec_j)), shape=(num_features, num_features))
    n_components, labels = connected_components(csgraph=mat_w, directed=False, return_labels=True)

    label_counts = Counter(labels).most_common()
    outlier_ratio = (num_features - label_counts[0][1]) / num_features

    return outlier_ratio