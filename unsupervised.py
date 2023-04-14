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