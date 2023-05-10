import os
HOME = os.path.expanduser('~')
import numpy as np
import random
random.seed(10)
from random import sample
from scipy.spatial.distance import pdist, cdist, squareform
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from collections import Counter
import igraph as ig
import leidenalg as la
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import PrecisionRecallDisplay, precision_recall_curve, auc
# import matplotlib
# matplotlib.use('tkagg')
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

def debug_features(mat_train, vec_labels, ref_motif):
    import matplotlib
    # matplotlib.use('tkagg')
    # import matplotlib.pyplot as plt
    import pickle
    from joblib import load

    classifier = 'logistic_regression'
    classifier_model_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/MAFIA_classifiers/A_m6A_multiple_noNorm_allReads'
    NUM_TOP_FEATURES = 10
    with open(os.path.join(classifier_model_dir, 'top_{}_model_weights.pkl'.format(NUM_TOP_FEATURES)), 'rb') as p_in:
        dict_top_feature_indices = pickle.load(p_in)
    top_feature_indices = dict_top_feature_indices[ref_motif]
    clf = load(os.path.join(classifier_model_dir, '{}_{}.joblib'.format(classifier, ref_motif)))
    probs = clf.predict_proba(mat_train)[:, 1]

    unm_features = mat_train[vec_labels==0]
    mod_features = mat_train[vec_labels==1]
    unm_probs = probs[vec_labels==0]
    mod_probs = probs[vec_labels==1]

    num_rows = 5
    plt.figure(figsize=(16, 10))
    for i in range(num_rows):
        plt.subplot(num_rows, 2, 2*i+1)
        plt.plot(unm_features[i], label='p={:.2f}'.format(unm_probs[i]))
        plt.plot(top_feature_indices, unm_features[i, top_feature_indices], 'ro')
        plt.xlim([0, len(unm_features[i])])
        plt.ylim([-0.5, 4])
        plt.legend(loc='upper left', fontsize=10)
        if i==0:
            plt.title('A', fontsize=15)
        if i==num_rows-1:
            plt.xlabel('Features', fontsize=15)
        else:
            plt.xticks([])

        plt.subplot(num_rows, 2, 2*i+2)
        plt.plot(mod_features[i], label='p={:.2f}'.format(mod_probs[i]))
        plt.plot(top_feature_indices, mod_features[i, top_feature_indices], 'ro')
        plt.xlim([0, len(mod_features[i])])
        plt.ylim([-0.5, 4])
        plt.legend(loc='upper left', fontsize=10)
        if i==0:
            plt.title('m6A', fontsize=15)
        if i==num_rows-1:
            plt.xlabel('Features', fontsize=15)
        else:
            plt.xticks([])

    plt.suptitle('{}, top {} features'.format(ref_motif, NUM_TOP_FEATURES), fontsize=20)
    plt.savefig(os.path.join(classifier_model_dir, '{}_top_{}_features.png'.format(ref_motif, NUM_TOP_FEATURES)), bbox_inches='tight')
    plt.close('all')

def train_binary_classifier(unm_dict, mod_dict, classifier, scaler=None, frac_test_split=0.25, debug_img_path=None, fig_title=None):
    import matplotlib
    # matplotlib.use('tkagg')
    # import matplotlib.pyplot as plt

    max_num_features = min(len(unm_dict), len(mod_dict))
    sample_unm_dict = {k: unm_dict[k] for k in sample(unm_dict.keys(), max_num_features)}
    sample_mod_dict = {k: mod_dict[k] for k in sample(mod_dict.keys(), max_num_features)}

    labels = np.array([0 for ii in range(len(sample_unm_dict))] + [1 for ii in range(len(sample_mod_dict))])
    unm_motifs = [v[0] for k, v in sample_unm_dict.items()]
    mod_motifs = [v[0] for k, v in sample_mod_dict.items()]
    unm_features = [v[1] for k, v in sample_unm_dict.items()]
    mod_features = [v[1] for k, v in sample_mod_dict.items()]
    all_features = np.array(unm_features + mod_features)

    X_train, X_test, y_train, y_test = train_test_split(all_features, labels, test_size=frac_test_split)
    # debug_features(X_train, y_train, ref_motif)

    if classifier=='svm':
        clf = SVC(gamma='auto', random_state=0, max_iter=1000)
    elif classifier=='logistic_regression':
        clf = LogisticRegression(random_state=0, max_iter=1000)

    if scaler=='MaxAbs':
        binary_model = make_pipeline(MaxAbsScaler(), clf)
    elif scaler=='Standard':
        binary_model = make_pipeline(StandardScaler(), clf)
    else:
        binary_model = clf

    binary_model = binary_model.fit(X_train, y_train)
    # predictions = binary_model.predict(X_test)
    # accuracy = binary_model.score(X_test, y_test)

    y_score = binary_model.decision_function(X_test)
    precision, recall, thresholds = precision_recall_curve(y_test, y_score)
    score_auc = auc(recall, precision)
    opt_ind = np.where(precision >= 0.9)[0][0]
    opt_thresh = thresholds[opt_ind-1]

    ### debug ############################################################
    if debug_img_path is not None:
        debug_img_dir = os.path.dirname(debug_img_path)
        if not os.path.exists(debug_img_dir):
            os.makedirs(debug_img_dir, exist_ok=True)
        opt_recall = recall[opt_ind]
        opt_predictions = np.int32(y_score>=opt_thresh)
        opt_accuracy = np.mean(opt_predictions==y_test)

        plt.figure(figsize=(5, 5))
        plt.plot(recall, precision, '-')
        # plt.axvline(opt_recall, c='g')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0, 1.05])
        plt.title('{}\n{} unm, {} mod, test split {:.2f}\n AUC = {:.2f}'.format(fig_title, len(unm_features), len(mod_features), frac_test_split, score_auc))

        plt.savefig(debug_img_path, bbox_inches='tight')
        plt.close('all')

        ### output precision, recall, thresholds ###
        thresh_precision_recall = np.vstack((np.concatenate([thresholds, [100.0]]), precision, recall)).T
        np.savetxt(debug_img_path.replace('.png', '.txt'), thresh_precision_recall, fmt='%.2f', delimiter='\t', header='threshold\tprecision\trecall')
        np.savetxt(debug_img_path.replace('.png', '_numFeatures.txt'), [len(y_train), len(y_test)], delimiter='\t', header='Train features\tTest features')
    ######################################################################

    return score_auc, binary_model, opt_thresh

def train_svm_ivt_wt(ivt_dict, wt_dict, wanted_motif, site, debug_img_dir=None):
    frac_test_split = 0.25
    max_num_features = min(len(ivt_dict), len(wt_dict))
    sample_ivt_dict = {k: ivt_dict[k] for k in sample(ivt_dict.keys(), max_num_features)}
    sample_wt_dict = {k: wt_dict[k] for k in sample(wt_dict.keys(), max_num_features)}

    labels = np.array([0 for ii in range(len(sample_ivt_dict))] + [1 for ii in range(len(sample_wt_dict))])
    ivt_motifs = [v[0] for k, v in sample_ivt_dict.items()]
    wt_motifs = [v[0] for k, v in sample_wt_dict.items()]
    ivt_features = [v[1] for k, v in sample_ivt_dict.items()]
    wt_features = [v[1] for k, v in sample_wt_dict.items()]
    all_features = np.array(ivt_features + wt_features)

    if Counter(ivt_motifs).most_common(1)[0][0] != wanted_motif:
        print('IVT motif {} doesn\'t match reference one {}'.format(Counter(ivt_motifs).most_common(1)[0][0], wanted_motif))
        return -1

    X_train, X_test, y_train, y_test = train_test_split(all_features, labels, test_size=frac_test_split)
    # clf = make_pipeline(StandardScaler(), SVC(gamma='auto'))
    clf = SVC(gamma='auto')
    clf = clf.fit(X_train, y_train)
    # predictions = clf.predict(X_test)
    # accuracy = clf.score(X_test, y_test)

    y_score = clf.decision_function(X_test)
    precision, recall, thresholds = precision_recall_curve(y_test, y_score)
    score_auc = auc(recall, precision)
    opt_ind = np.where(precision >= 0.8)[0][0]
    opt_thresh = thresholds[opt_ind-1]

    ### debug ############################################################
    if debug_img_dir is not None:
        if not os.path.exists(debug_img_dir):
            os.makedirs(debug_img_dir, exist_ok=True)
        opt_recall = recall[opt_ind]
        opt_predictions = np.int32(y_score>=opt_thresh)
        opt_accuracy = np.mean(opt_predictions==y_test)

        plt.figure(figsize=(5, 5))
        plt.plot(recall, precision)
        plt.axvline(opt_recall, c='g')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0, 1.05])
        plt.title('{} {}\n{} WT, {} IVT, AUC = {:.2f}'.format(site['mod'], site['stop'], len(wt_features), len(ivt_features), score_auc))

        plt.savefig(os.path.join(debug_img_dir, 'svm_auc_{}_{}.png'.format(site['mod'], site['stop'])), bbox_inches='tight')
        plt.close('all')
        # plt.show()
        # display = PrecisionRecallDisplay.from_estimator(clf, X_test, y_test)
    ######################################################################

    return score_auc, clf, opt_thresh

def get_mod_ratio_with_binary_classifier(dict_motif_feature, clf, mod_thresh=0.5, output_mod_probs=False):
    test_motifs = [v[0] for k, v in dict_motif_feature.items()]
    test_features = [v[1] for k, v in dict_motif_feature.items()]

    mod_probs = clf.predict_proba(test_features)[:, 1]
    read_predMotif_modProbs = {read_id: (pred_motif, mod_prob) for read_id, pred_motif, mod_prob in zip(dict_motif_feature.keys(), test_motifs, mod_probs)}
    predictions = np.int32(mod_probs > mod_thresh)
    avg_mod_ratio = np.mean(predictions)

    if output_mod_probs:
        return avg_mod_ratio, read_predMotif_modProbs
    else:
        return avg_mod_ratio

def get_mod_probs(dict_motif_feature, clf):
    test_motifs = [v[0] for k, v in dict_motif_feature.items()]
    test_features = [v[1] for k, v in dict_motif_feature.items()]
    mod_probs = clf.predict_proba(test_features)[:, 1]
    read_mod_probs = {k: v for k, v in zip(dict_motif_feature.keys(), mod_probs)}
    return read_mod_probs

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