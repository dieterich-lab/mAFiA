import os
import numpy as np
from sklearn.model_selection import train_test_split
from umap import UMAP
from collections import Counter

def plot_umap(in_embedding, labels, title, img_out=None):
    # import matplotlib
    # matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    labels = np.array(labels)

    plt.figure(figsize=(8, 8))
    for this_label in np.unique(labels):
        sub_embedding = in_embedding[labels==this_label]
        plt.scatter(sub_embedding[:, 0], sub_embedding[:, 1], label=this_label)
    plt.legend(fontsize=15)
    plt.xlabel('UMAP 1', fontsize=15)
    plt.ylabel('UMAP 2', fontsize=15)
    plt.title(title, fontsize=20)
    if img_out is not None:
        if not os.path.exists(img_out):
            os.makedirs(img_out, exist_ok=True)
        plt.savefig(os.path.join(img_out, 'umap_{}.png'.format('_'.join(title.split(' ')))), bbox_inches='tight')
        plt.close('all')

def train_cluster(unm_dict, mod_dict, ref_motif, scaler=None, debug_img_dir=None):
    print('Now clustering features...', flush=True)
    labels = np.array([0 for ii in range(len(unm_dict))] + [1 for ii in range(len(mod_dict))])
    unm_motifs = [v[0] for k, v in unm_dict.items()]
    mod_motifs = [v[0] for k, v in mod_dict.items()]
    unm_features = [v[1] for k, v in unm_dict.items()]
    mod_features = [v[1] for k, v in mod_dict.items()]
    all_features = np.array(unm_features + mod_features)

    # frac_test_split = 0.25
    # X_train, X_test, y_train, y_test = train_test_split(all_features, labels, test_size=frac_test_split)
    # # debug_features(X_train, y_train, ref_motif)

    if scaler=='MaxAbs':
        scaled_features = all_features / np.max(all_features, axis=1)
    else:
        scaled_features = all_features

    ######################################################################################################
    ### umap ###
    reducer = UMAP(random_state=0)
    embedding = reducer.fit_transform(scaled_features)
    label_names = ['unm' if this_label==0 else 'mod' for this_label in labels]
    plot_umap(embedding, label_names, '{} cluster'.format(ref_motif), debug_img_dir)

    ######################################################################################################
    # if classifier=='svm':
    #     clf = SVC(gamma='auto', random_state=0, max_iter=1000)
    # elif classifier=='logistic_regression':
    #     clf = LogisticRegression(random_state=0, max_iter=1000)
    #
    # if scaler is None:
    #     binary_model = clf
    # elif scaler=='MaxAbs':
    #     binary_model = make_pipeline(MaxAbsScaler(), clf)
    # elif scaler=='Standard':
    #     binary_model = make_pipeline(StandardScaler(), clf)
    #
    # binary_model = binary_model.fit(X_train, y_train)
    # # predictions = binary_model.predict(X_test)
    # # accuracy = binary_model.score(X_test, y_test)
    #
    # y_score = binary_model.decision_function(X_test)
    # precision, recall, thresholds = precision_recall_curve(y_test, y_score)
    # score_auc = auc(recall, precision)
    # opt_ind = np.where(precision >= 0.9)[0][0]
    # opt_thresh = thresholds[opt_ind-1]
    #
    # ### debug ############################################################
    # if debug_img_dir is not None:
    #     if not os.path.exists(debug_img_dir):
    #         os.makedirs(debug_img_dir, exist_ok=True)
    #     opt_recall = recall[opt_ind]
    #     opt_predictions = np.int32(y_score>=opt_thresh)
    #     opt_accuracy = np.mean(opt_predictions==y_test)
    #
    #     plt.figure(figsize=(5, 5))
    #     plt.plot(recall, precision)
    #     # plt.axvline(opt_recall, c='g')
    #     plt.xlabel('Recall')
    #     plt.ylabel('Precision')
    #     plt.ylim([0, 1.05])
    #     plt.title('{}\n{} unm, {} mod, AUC = {:.2f}'.format(ref_motif, len(unm_features), len(mod_features), score_auc))
    #
    #     plt.savefig(os.path.join(debug_img_dir, 'log_reg_auc_{}.png'.format(ref_motif)), bbox_inches='tight')
    #     plt.close('all')
    # ######################################################################
    #
    # return score_auc, binary_model, opt_thresh

    return 1