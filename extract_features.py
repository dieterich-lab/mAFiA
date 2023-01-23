import os, sys
HOME = os.path.expanduser('~')
from utils import segment
from tqdm import tqdm
import torch
import numpy as np
from Bio.Seq import Seq
from utils import get_norm_signal_from_read_id
from models import objectview, rodan
from fast_ctc_decode import beam_search

model_path = os.path.join(HOME, 'pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch')
torchdict = torch.load(model_path, map_location="cpu")
origconfig = torchdict["config"]
config = objectview(origconfig)

def convert_statedict(state_dict):
    from collections import OrderedDict
    new_checkpoint = OrderedDict()
    for k, v in state_dict.items():
        name = k[7:] # remove module.
        new_checkpoint[name] = v
    return new_checkpoint

activation = {}
def get_activation(name):
    def hook(model, input, output):
        activation[name] = output.detach()
    return hook

def load_model(modelfile, config):
    if modelfile == None:
        sys.stderr.write("No model file specified!")
        sys.exit(1)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # print("Using device:", device)
    model = rodan(config=config).to(device)
    state_dict = torch.load(modelfile, map_location=device)["state_dict"]
    if "state_dict" in state_dict:
        model.load_state_dict(convert_statedict(state_dict["state_dict"]))
    else:
        model.load_state_dict(torch.load(modelfile, map_location=device)["state_dict"])
    model.convlayers.conv21.register_forward_hook(get_activation('conv21'))

    # print(model)
    model.eval()
    torch.set_grad_enabled(False)

    return model, device

vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }
alphabet = "".join(["N"] + list(vocab.values()))

def ctcdecoder(logits, label, blank=False, beam_size=5, alphabet=alphabet, pre=None):
    ret = np.zeros((label.shape[0], label.shape[1]+50))
    retstr = []
    for i in range(logits.shape[0]):
        if pre is not None:
            beamcur = beam_search(torch.softmax(torch.tensor(pre[:,i,:]), dim=-1).cpu().detach().numpy(), alphabet=alphabet, beam_size=beam_size)[0]
        prev = None
        cur = []
        pos = 0
        for j in range(logits.shape[1]):
            if not blank:
                if logits[i,j] != prev:
                    prev = logits[i,j]
                    try:
                        if prev != 0:
                            ret[i, pos] = prev
                            pos+=1
                            cur.append(vocab[prev])
                    except:
                        sys.stderr.write("ctcdecoder: fail on i:", i, "pos:", pos)
            else:
                if logits[i,j] == 0: break
                ret[i, pos] = logits[i,j] # is this right?
                cur.append(vocab[logits[i,pos]])
                pos+=1
        if pre is not None:
            retstr.append(beamcur)
        else:
            retstr.append("".join(cur))
    return ret, retstr

def extract_features_from_signal(model, device, signal, pos, bam_motif):
    chunks = segment(signal, config.seqlen)
    event = torch.unsqueeze(torch.FloatTensor(chunks), 1).to(device, non_blocking=True)
    out = model.forward(event)

    pred_labels = []
    features = []
    for i in range(out.shape[1]):
        this_pred_label, pred_locs = beam_search(torch.softmax(out[:, i, :], dim=-1).cpu().detach().numpy(), alphabet=alphabet)
        pred_labels.append(this_pred_label)
        this_feature = activation['conv21'].detach().cpu().numpy()[i, :, pred_locs]
        features.append(this_feature)

    pred_label = ''.join(pred_labels)[::-1]
    all_features = np.vstack(features)[::-1]

    pred_motif = pred_label[pos-2:pos+3]
    # print('Predicted motif {}, aligned {}, matching {}'.format(pred_motif, bam_motif, pred_motif == bam_motif))
    if pred_motif!=bam_motif:
        print('\n!!!!!!!!!!!!!!!!!!!Error: Predicted motif =/= aligned!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        return None

    return all_features[pos], pred_motif

def collect_features_from_aligned_site(alignment, index_read_ids, chr, site):
    fixed_model, fixed_device = load_model(model_path, config)

    site_motif_features = {}
    for pileupcolumn in alignment.pileup(chr, site, site+1, truncate=True, min_base_quality=0):
        if pileupcolumn.pos == site:
            coverage = pileupcolumn.get_num_aligned()
            if coverage>100:
                valid_counts = 0
                for pileupread in tqdm(pileupcolumn.pileups):
                    query_name = pileupread.alignment.query_name
                    query_position = pileupread.query_position_or_next
                    # query_position = pileupread.query_position
                    flag = pileupread.alignment.flag

                    # if query_position and (flag==0) and (query_name in index_read_ids.keys()) and (pileupread.alignment.query_sequence[query_position] == 'A'):
                    if query_position and (flag==0) and (query_name in index_read_ids.keys()):
                    # if query_position and (flag==0) and (query_name in id_signal.keys()):
                        valid_counts += 1
                        query_motif = pileupread.alignment.query_sequence[query_position-2:query_position+3]
                        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
                        # this_read_signal = id_signal[query_name]
                        this_read_features, this_read_motif = extract_features_from_signal(fixed_model, fixed_device, this_read_signal, query_position, query_motif)
                        if this_read_features is not None:
                            site_motif_features[query_name] = (this_read_motif, this_read_features)
    return site_motif_features