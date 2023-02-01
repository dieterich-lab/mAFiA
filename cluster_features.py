import os
HOME = os.path.expanduser('~')
import numpy as np
from scipy.spatial.distance import pdist, squareform
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
    # label_counts = Counter(labels).most_common()
    # outlier_ratio = (dim - label_counts[0][1]) / dim

    return labels

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
    # label_counts = Counter(partition.membership).most_common()
    # outlier_ratio = (dim - label_counts[0][1]) / dim

    return partition.membership

def get_outlier_ratio_from_features(ivt_dict, wt_dict, wanted_motif, perc_thresh=0.8):
    labels = ['ivt' for ii in range(len(ivt_dict))] + ['wt' for ii in range(len(wt_dict))]
    ivt_motifs = [v[0] for k, v in ivt_dict.items()]
    wt_motifs = [v[0] for k, v in wt_dict.items()]

    if Counter(ivt_motifs).most_common(1)[0][0] != wanted_motif:
        print('IVT motif {} doesn\'t match reference one {}'.format(Counter(ivt_motifs).most_common(1)[0][0], wanted_motif))
        return -1

    all_dicts = dict(ivt_dict)
    all_dicts.update(wt_dict)
    all_ids = []
    all_motifs = []
    all_features = []
    for k, v in all_dicts.items():
        all_ids.append(k)
        all_motifs.append(v[0])
        all_features.append(v[1])

    ### check features ###
    # plt.figure(figsize=(16, 16))
    # for i in range(16):
    #     plt.subplot(4, 4, i+1)
    #     plt.plot(all_features[i])
    # plt.savefig(os.path.join(img_out, 'feature_vectors.png'), bbox_inches='tight')
    # plt.close('all')

    num_features = len(all_features)
    vec_w = 1.0 - pdist(np.vstack(all_features), metric='cosine')
    # mat_w = squareform(vec_w)

    # for perc_thresh in np.arange(0, 1, 0.01):
        # membership = cluster_by_louvain(vec_w, num_features)

    membership = cluster_by_connected_components(vec_w, num_features, perc_thresh)
    # if len(np.unique(ivt_membership))>30:
    #     break
    # motifs = np.array(ivt_motifs)[ivt_membership==0]
    # print(motifs, Counter(motifs).most_common())
    # if (len(Counter(motifs).most_common())==1) and Counter(motifs).most_common()[0][0]==wanted_motif:
    #     break

    # print(perc_thresh, np.unique(ivt_membership))
    ivt_membership = membership[np.array(labels)=='ivt']
    wt_membership = membership[np.array(labels)=='wt']
    outlier_ratio = 1.0 - sum(wt_membership==0) / len(wt_membership)
    # print('Outlier ratio = {:.2f}'.format(outlier_ratio))

    return outlier_ratio

    # plt.figure(figsize=(8, 8))
    # plt.hist(vec_w, bins=50)
    # plt.savefig(os.path.join(img_out, 'hist_cosine_similarities.png'), bbox_inches='tight')
    # plt.close('all')