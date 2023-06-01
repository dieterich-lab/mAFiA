import os
import numpy as np
from umap import UMAP
from collections import Counter
from scipy.spatial.distance import pdist, cdist, squareform
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import igraph as ig
import leidenalg as la
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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

def plot_umap(in_embedding, labels, title, img_out=None):
    labels = np.array(labels)

    plt.figure(figsize=(8, 8))
    for this_label in np.unique(labels):
        sub_embedding = in_embedding[labels==this_label]
        plt.scatter(sub_embedding[:, 0], sub_embedding[:, 1], label='{} ({})'.format(this_label, np.sum(labels==this_label)))
    plt.legend(fontsize=15)
    plt.xlabel('UMAP 1', fontsize=15)
    plt.ylabel('UMAP 2', fontsize=15)
    plt.title(title, fontsize=20)
    if img_out is not None:
        if not os.path.exists(img_out):
            os.makedirs(img_out, exist_ok=True)
        plt.savefig(os.path.join(img_out, 'umap_{}.png'.format(title.replace(' ', '_').replace('\n', '_'))), bbox_inches='tight')
        plt.close('all')

def train_cluster(unm_dict, mod_dict, site_name, scaler=None, debug_img_dir=None):
    print('Now clustering features...', flush=True)
    labels = np.array([0 for ii in range(len(unm_dict))] + [1 for ii in range(len(mod_dict))])
    unm_motifs = [v[0] for k, v in unm_dict.items()]
    mod_motifs = [v[0] for k, v in mod_dict.items()]
    unm_features = [v[1] for k, v in unm_dict.items()]
    mod_features = [v[1] for k, v in mod_dict.items()]
    all_features = np.array(unm_features + mod_features)

    if scaler=='MaxAbs':
        scaled_features = all_features / np.max(all_features, axis=1)
    else:
        scaled_features = all_features

    ### umap ###
    reducer = UMAP(random_state=0)
    embedding = reducer.fit_transform(scaled_features)
    label_names = ['IVT' if this_label==0 else 'WT' for this_label in labels]
    plot_umap(embedding, label_names, '{}'.format(site_name), debug_img_dir)

    return 1

def calculate_outlier_ratio_with_ivt_distance(ivt_dict, wt_dict, scaler=None, sigma_multiplier=1, debug_img_dir=None):
    print('Now clustering features...', flush=True)
    labels = np.array([0 for ii in range(len(ivt_dict))] + [1 for ii in range(len(wt_dict))])
    ivt_features = np.vstack([v[1] for k, v in ivt_dict.items()])
    wt_features = np.vstack([v[1] for k, v in wt_dict.items()])

    if scaler=='MaxAbs':
        print('Re-scaling features with {}'.format(scaler))
        ivt_features = ivt_features / np.max(ivt_features, axis=1, keepdims=True)
        wt_features = wt_features / np.max(wt_features, axis=1, keepdims=True)

    ivt_centroid = np.mean(ivt_features, axis=0)
    ivt_l1_dist = np.sum(np.abs(ivt_features - ivt_centroid), axis=1)
    ivt_l1_dist_mu = np.mean(ivt_l1_dist)
    ivt_l1_dist_sigma = np.std(ivt_l1_dist)

    wt_l1_dist = np.sum(np.abs(wt_features - ivt_centroid), axis=1)
    pred_mod_ratio = np.mean(np.abs(wt_l1_dist - ivt_l1_dist_mu) > (ivt_l1_dist_sigma * sigma_multiplier)) * 100.0

    return pred_mod_ratio

def get_outlier_ratio_from_features(ivt_dict, wt_dict, wanted_motif, perc_thresh=0.8):
    labels = ['ivt' for ii in range(len(ivt_dict))] + ['wt' for ii in range(len(wt_dict))]
    ivt_motifs = [v[0] for k, v in ivt_dict.items()]
    wt_motifs = [v[0] for k, v in wt_dict.items()]
    ivt_features = [v[1] for k, v in ivt_dict.items()]
    wt_features = [v[1] for k, v in wt_dict.items()]

    if Counter(ivt_motifs).most_common(1)[0][0] != wanted_motif:
        print('IVT motif {} doesn\'t match reference one {}'.format(Counter(ivt_motifs).most_common(1)[0][0], wanted_motif))
        return -1
    # if Counter(wt_motifs).most_common(1)[0][0] != wanted_motif:
    #     print('WT motif {} doesn\'t match reference one {}'.format(Counter(wt_motifs).most_common(1)[0][0], wanted_motif))
    #     return -1

    all_dicts = dict(ivt_dict)
    all_dicts.update(wt_dict)
    all_ids = []
    all_motifs = []
    all_features = []
    for k, v in all_dicts.items():
        all_ids.append(k)
        all_motifs.append(v[0])
        all_features.append(v[1])

    ### debug ###
    # import matplotlib
    # matplotlib.use('tkagg')
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(16, 16))
    # for i in range(8):
    #     plt.subplot(8, 2, 2*i+1)
    #     plt.plot(ivt_features[i])
    #     if i==0:
    #         plt.title('IVT')
    #     plt.xlim([-1, 769])
    #     if i<7:
    #         plt.xticks([])
    #     else:
    #         plt.xlabel('Feature')
    #     plt.subplot(8, 2, 2*i+2)
    #     plt.plot(wt_features[i])
    #     if i==0:
    #         plt.title('WT')
    #     plt.xlim([-1, 769])
    #     if i<7:
    #         plt.xticks([])
    #     else:
    #         plt.xlabel('Feature')

    # plt.savefig(os.path.join(img_out, 'feature_vectors.png'), bbox_inches='tight')
    # plt.close('all')
    ##################################

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
    ivt_largest_cluster = Counter(ivt_membership).most_common()[0][0]
    # if ivt_largest_cluster!=0:
    #     print('IVT largest cluster [{}] not 0!'.format(ivt_largest_cluster))
    #     return -1
    wt_membership = membership[np.array(labels)=='wt']
    outlier_ratio = 1.0 - sum(wt_membership==ivt_largest_cluster) / len(wt_membership)
    # print('Outlier ratio = {:.2f}'.format(outlier_ratio))

    return outlier_ratio

    # plt.figure(figsize=(8, 8))
    # plt.hist(vec_w, bins=50)
    # plt.savefig(os.path.join(img_out, 'hist_cosine_similarities.png'), bbox_inches='tight')
    # plt.close('all')

def calculate_self_intersection(indices):
    num_indices = len(indices)
    num_intersections = []
    for i in range(num_indices):
        for j in range(i+1, num_indices):
            num_intersections.append(len(indices[i].intersection(indices[j])))
    return num_intersections

def calculate_cross_intersection(indices_base, indices_compare):
    num_intersections = []
    for i in range(len(indices_base)):
        for j in range(len(indices_compare)):
            num_intersections.append(len(indices_base[i].intersection(indices_compare[j])))
    return num_intersections

def get_outlier_ratio_from_features_v2(ivt_dict, wt_dict, wanted_motif, sigma_dev = 1.0):
    ivt_motifs = [v[0] for k, v in ivt_dict.items()]
    wt_motifs = [v[0] for k, v in wt_dict.items()]
    if Counter(ivt_motifs).most_common(1)[0][0] != wanted_motif:
        print('IVT motif {} doesn\'t match reference one {}'.format(Counter(ivt_motifs).most_common(1)[0][0], wanted_motif))
        return -1

    ivt_features = np.vstack([v[1] for k, v in ivt_dict.items()])
    # ivt_features -= np.mean(ivt_features, axis=1)[:, np.newaxis]
    # ivt_features /= np.std(ivt_features, axis=1)[:, np.newaxis]
    wt_features = np.vstack([v[1] for k, v in wt_dict.items()])
    # wt_features -= np.mean(wt_features, axis=1)[:, np.newaxis]
    # wt_features /= np.std(wt_features, axis=1)[:, np.newaxis]

    c_mat_ivt = np.mean(ivt_features[:, np.newaxis, :] * ivt_features[np.newaxis, :, ], axis=-1)
    i_ind, j_ind = np.triu_indices_from(c_mat_ivt, k=1)
    c_vec_ivt = c_mat_ivt[i_ind, j_ind]
    c_mat_wt_vs_ivt = np.mean(wt_features[:, np.newaxis, :] * ivt_features[np.newaxis, :, ], axis=-1)

    dev_vec = (np.mean(c_mat_wt_vs_ivt, axis=1) - np.mean(c_vec_ivt)) / np.std(c_vec_ivt)
    num_outliers = np.sum(dev_vec<=-sigma_dev)
    return num_outliers / len(wt_motifs)

def get_outlier_ratio_from_features_v3(ivt_dict, wt_dict, wanted_motif, sigma_dev=1.0):
    ivt_motifs = [v[0] for k, v in ivt_dict.items()]
    wt_motifs = [v[0] for k, v in wt_dict.items()]
    if Counter(ivt_motifs).most_common(1)[0][0] != wanted_motif:
        print('IVT motif {} doesn\'t match reference one {}'.format(Counter(ivt_motifs).most_common(1)[0][0], wanted_motif))
        return -1

    ivt_features = np.vstack([v[1] for k, v in ivt_dict.items()])
    wt_features = np.vstack([v[1] for k, v in wt_dict.items()])

    c_vec_ivt = 1.0 - pdist(ivt_features, metric='cosine')
    c_vec_ivt_wt = 1.0 - pdist(np.vstack([ivt_features, wt_features]), metric='cosine')
    c_mat_ivt_wt = squareform(c_vec_ivt_wt)
    c_mat_wt_vs_ivt = c_mat_ivt_wt[len(ivt_motifs):, :][:, :len(ivt_motifs)]

    dev_vec = (np.mean(c_mat_wt_vs_ivt, axis=1) - np.mean(c_vec_ivt)) / np.std(c_vec_ivt)
    num_outliers = np.sum(dev_vec<=-sigma_dev)
    return num_outliers / len(wt_motifs)