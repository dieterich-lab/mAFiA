import os, sys
HOME = os.path.expanduser('~')
sys.path.append(os.path.join(HOME, 'git/mAFiA_dev'))
import pickle
from mAFiA.feature_classifiers import MotifClassifier

classifier_type = 'logistic_regression'
scaler = 'MaxAbs'

prj_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A'

read_id_paths_unm = [
    os.path.join(prj_dir, 'oligo/read_ids_ISA-WUE_A.txt'),
    os.path.join(prj_dir, 'oligo/read_ids_ISA_mix1-4_A.txt'),
    os.path.join(prj_dir, 'oligo/read_ids_ISA_mix17-19_A.txt')
]

read_id_paths_mod = [
    os.path.join(prj_dir, 'oligo/read_ids_ISA-WUE_m6A.txt'),
    os.path.join(prj_dir, 'oligo/read_ids_ISA_mix1-4_m6A.txt'),
    os.path.join(prj_dir, 'oligo/read_ids_ISA_mix20-22_m6A.txt')
]

### AAACA ###
# old_clf_paths = [
#     os.path.join(prj_dir, 'MAFIA_classifiers/_DRACH/AAACA.pkl'),
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix17_mix20/AAACA.pkl'),
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix19_mix22/AAACA.pkl')
# ]

### GAACT ###
old_clf_paths = [
    os.path.join(prj_dir, 'MAFIA_classifiers/ISA-WUE/GAACT.pkl'),
    os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix17_mix20/GAACT.pkl')
]

### TAACC ###
# old_clf_paths = [
#     os.path.join(prj_dir, 'MAFIA_classifiers/_DRACH/TAACC.pkl'),
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix18_mix21/TAACC.pkl')
# ]

### TAACT ###
# old_clf_paths = [
#     os.path.join(prj_dir, 'MAFIA_classifiers/_DRACH/TAACT.pkl'),
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix17_mix20/TAACT.pkl'),
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix18_mix21/TAACT.pkl')
# ]

### TGACT ###
# old_clf_paths = [
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA-WUE/TGACT.pkl'),
#     os.path.join(prj_dir, 'MAFIA_classifiers/ISA_mix19_mix22/TGACT.pkl')
# ]

out_dir = os.path.join(prj_dir, 'MAFIA_classifiers/DRACH5_retrain')
os.makedirs(out_dir, exist_ok=True)

read_ids_unm = []
for this_path in read_id_paths_unm:
    with open(this_path, 'r') as h:
        read_ids_unm.extend([l.rstrip('\n') for l in h.readlines()])

read_ids_mod = []
for this_path in read_id_paths_mod:
    with open(this_path, 'r') as h:
        read_ids_mod.extend([l.rstrip('\n') for l in h.readlines()])

dict_labels = {}
for this_id in read_ids_unm:
    dict_labels[this_id] = 0
for this_id in read_ids_mod:
    dict_labels[this_id] = 1

motif = [os.path.basename(path).split('.')[0] for path in old_clf_paths]
if len(set(motif))>1:
    print('Warning: non-unique motif!')
    exit(0)
else:
    motif = motif[0]

all_nts = []
for this_path in old_clf_paths:
    with open(this_path, 'rb') as h:
        this_clf = pickle.load(h)
    all_nts.append(this_clf.train_nts)
    all_nts.append(this_clf.test_nts)
all_nts = [this_nt for ll in all_nts for this_nt in ll]
all_labels = [dict_labels.get(this_nt.read_id, -1) for this_nt in all_nts]

unm_nts = []
mod_nts = []
for this_nt, this_label in zip(all_nts, all_labels):
    if this_label==0:
        unm_nts.append(this_nt)
    elif this_label==1:
        mod_nts.append(this_nt)

########################################################################################################################
### Train new classifier ###############################################################################################
########################################################################################################################
new_classifier = MotifClassifier(motif=motif, classifier_type=classifier_type, scaler=scaler)
new_classifier.train(unm_nts, mod_nts, frac_test_split=0.2)
new_classifier.save(
    out_model_path=os.path.join(out_dir, '{}.pkl'.format(motif)),
    draw_prc=True
)
