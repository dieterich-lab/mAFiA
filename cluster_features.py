import os
HOME = os.path.expanduser('~')
import numpy as np
from scipy.spatial.distance import pdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from collections import Counter
import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt

img_out = os.path.join(HOME, 'img_out')

def cluster_by_connected_components(vec_s, dim, threshold=0.5):
    vec_i, vec_j = np.triu_indices(dim, k=1)
    mask = (vec_s >= threshold)
    sel_vec_i = vec_i[mask]
    sel_vec_j = vec_j[mask]
    sel_vec_w = vec_s[mask]
    mat_w = coo_matrix((sel_vec_w, (sel_vec_i, sel_vec_j)), shape=(dim, dim))
    n_components, labels = connected_components(csgraph=mat_w, directed=False, return_labels=True)
    label_counts = Counter(labels).most_common()
    outlier_ratio = (dim - label_counts[0][1]) / dim

    return outlier_ratio

def cluster_by_louvain(vec_s, dim):
    # vec_s_zeroed = vec_s - np.min(vec_s)
    vec_s_zeroed = vec_s.copy()
    vec_s_zeroed[vec_s_zeroed<0] = 0
    vec_i, vec_j = np.triu_indices(dim, k=1)
    g = ig.Graph()
    g.add_vertices(dim)
    g.add_edges(zip(vec_i, vec_j))
    g.es['weight'] = vec_s_zeroed
    partition = la.ModularityVertexPartition(g, weights='weight')
    label_counts = Counter(partition.membership).most_common()
    outlier_ratio = (dim - label_counts[0][1]) / dim

    return outlier_ratio

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

    # calc_ratio = cluster_by_connected_components(vec_w, num_features)
    calc_ratio = cluster_by_louvain(vec_w, num_features)

    return calc_ratio

    # plt.figure(figsize=(8, 8))
    # plt.hist(vec_w, bins=50)
    # plt.savefig(os.path.join(img_out, 'hist_cosine_similarities.png'), bbox_inches='tight')
    # plt.close('all')