import os, sys
HOME = os.path.expanduser('~')
from utils import segment
from tqdm import tqdm
import torch
import numpy as np
from Bio.Seq import Seq
from utils import get_norm_signal_from_read_id
from models import objectview, rodan
from fast_ctc_decode import beam_search, viterbi_search
from class_defs import nucleotide, aligned_read
from time import time
import random
random.seed(10)
from random import sample

CTC_MODE='viterbi'
if CTC_MODE=='beam':
    ctc_decoder = beam_search
elif CTC_MODE=='viterbi':
    ctc_decoder = viterbi_search

def convert_statedict(state_dict):
    from collections import OrderedDict
    new_checkpoint = OrderedDict()
    for k, v in state_dict.items():
        name = k[7:] # remove module.
        new_checkpoint[name] = v
    return new_checkpoint

def get_freer_device():
    if torch.cuda.is_available():
        os.system('nvidia-smi -q -d Memory |grep -A5 GPU|grep Free >tmp')
        memory_available = [int(x.split()[2]) for x in open('tmp', 'r').readlines()]
        device_num = np.argmax(memory_available)
        freer_device = torch.device("cuda:{}".format(device_num))
        free_mem = memory_available[device_num]
    else:
        freer_device = torch.device("cpu")
        free_mem = -1
    return freer_device, free_mem

activation = {}
def get_activation(name):
    def hook(model, input, output):
        activation[name] = output.detach()
    return hook

def load_model(modelfile, config, ext_layer):
    if modelfile == None:
        sys.stderr.write("No model file specified!")
        sys.exit(1)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)
    # device, mem = get_freer_device()
    # print("Using device {} with free memory {}MB".format(device, mem), flush=True)
    model = rodan(config=config).to(device)
    state_dict = torch.load(modelfile, map_location=device)["state_dict"]
    if "state_dict" in state_dict:
        model.load_state_dict(convert_statedict(state_dict["state_dict"]))
    else:
        model.load_state_dict(torch.load(modelfile, map_location=device)["state_dict"])
    # model.convlayers.conv21.register_forward_hook(get_activation('conv21'))
    # model.convlayers.conv20.register_forward_hook(get_activation('conv20'))
    model.get_submodule(ext_layer).register_forward_hook(get_activation(ext_layer))

    # print(model)
    model.eval()
    torch.set_grad_enabled(False)

    return model, device

vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }
alphabet = "".join(["N"] + list(vocab.values()))
alphabet_to_num = {v: k for k, v in enumerate(list(alphabet))}

def ctcdecoder(logits, label, blank=False, beam_size=5, alphabet=alphabet, pre=None):
    ret = np.zeros((label.shape[0], label.shape[1]+50))
    retstr = []
    for i in range(logits.shape[0]):
        if pre is not None:
            beamcur = ctc_decoder(torch.softmax(torch.tensor(pre[:,i,:]), dim=-1).cpu().detach().numpy(), alphabet=alphabet, beam_size=beam_size)[0]
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

def get_features_from_signal(model, device, config, signal, ext_layer, feature_width=0):
    chunks = segment(signal, config.seqlen)
    event = torch.unsqueeze(torch.FloatTensor(chunks), 1).to(device, non_blocking=True)
    event_size = event.shape[0]
    batch_size = config.batchsize
    if event_size<=batch_size:
        out = model.forward(event)
        layer_activation = activation[ext_layer].detach().cpu().numpy()
    else:
        break_pts = np.arange(0, event_size, batch_size)
        start_stop_pts = [(start, stop) for start, stop in zip(break_pts, list(break_pts[1:]) + [event_size])]
        batch_out = []
        batch_layer_activation = []
        for (start, stop) in start_stop_pts:
            batch_out.append(model.forward(event[start:stop]))
            batch_layer_activation.append(activation[ext_layer].detach().cpu().numpy())
        out = torch.concat(batch_out, 1)
        layer_activation = np.concatenate(batch_layer_activation, axis=0)
    num_locs = layer_activation.shape[-1]
    pred_labels = []
    features = []
    for i in range(out.shape[1]):
        probs = torch.softmax(out[:, i, :], dim=-1).cpu().detach().numpy()
        this_pred_label, pred_locs = ctc_decoder(probs, alphabet=alphabet)

        pred_labels.append(this_pred_label)

        pred_locs_corrected = []
        for loc_ind in range(len(this_pred_label)):
            start_loc = pred_locs[loc_ind]
            if loc_ind<(len(this_pred_label)-1):
                end_loc = pred_locs[loc_ind+1]
            else:
                end_loc = probs.shape[0]
            pred_locs_corrected.append(start_loc + np.argmax(probs[start_loc:end_loc, alphabet_to_num[this_pred_label[loc_ind]]]))

        ### debug feature loc ###
        # import matplotlib
        # matplotlib.use('tkagg')
        # import matplotlib.pyplot as plt
        # alphabet_colors = {a: c for (a, c) in zip(alphabet, ['black', 'blue', 'red', 'green', 'purple'])}
        #
        # plt.figure(figsize=(20, 10))
        # plt.subplot(2, 1, 1)
        # for i in range(probs.shape[1]):
        #     plt.plot(probs[:, i], label=alphabet[i], color=alphabet_colors[alphabet[i]])
        #     plt.legend(loc='upper left')
        # for i, loc in enumerate(pred_locs):
        #     plt.axvline(loc, color=alphabet_colors[this_pred_label[i]])
        # plt.xticks(pred_locs, this_pred_label)
        # plt.title('Original')
        # plt.subplot(2, 1, 2)
        # for i in range(probs.shape[1]):
        #     plt.plot(probs[:, i], label=alphabet[i], color=alphabet_colors[alphabet[i]])
        #     plt.legend(loc='upper left')
        # for i, loc in enumerate(pred_locs_corrected):
        #     plt.axvline(loc, color=alphabet_colors[this_pred_label[i]])
        # plt.xticks(pred_locs_corrected, this_pred_label)
        # plt.title('Corrected')

        #########################

        if feature_width==0:
            this_feature = layer_activation[i, :, pred_locs_corrected]
        else:
            this_feature = []
            for loc_shift in range(-feature_width, feature_width+1):
                shifted_locs = [max(min(x+loc_shift, num_locs-1), 0) for x in pred_locs_corrected]
                this_feature.append(layer_activation[i, :, shifted_locs])
            this_feature = np.hstack(this_feature)
            # this_feature = np.mean(np.stack(this_feature, axis=-1), axis=-1)
        features.append(this_feature)
    pred_label = ''.join(pred_labels)[::-1]
    features = np.vstack(features)[::-1]
    return features, pred_label

def get_features_from_collection_of_signals(model, device, config, in_index_read_ids, max_num_reads, layer, feature_width=0):
    if max_num_reads>0:
        index_read_ids = {id: in_index_read_ids[id] for id in sample(list(in_index_read_ids.keys()), min(len(in_index_read_ids.keys()), max_num_reads))}
    else:
        index_read_ids = in_index_read_ids

    print('Extracting features from {}, width {}...'.format(layer, feature_width))
    id_predStr_feature = {}
    for query_name in tqdm(index_read_ids.keys()):
        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
        this_read_features, this_read_predStr = get_features_from_signal(model, device, config, this_read_signal, layer, feature_width)
        id_predStr_feature[query_name] = (this_read_predStr, this_read_features)
    return id_predStr_feature

########################################################################################################################
### oligo ##############################################################################################################
########################################################################################################################

def get_single_motif_nucleotides(motif_ind, reference, bam_in, predStr_features, block_size, block_center, enforce_ref_5mer=False):
    relevant_contigs = [k for k in reference.keys() if ((motif_ind in k.split('_')[1]) and (k in bam_in.references))]
    motif_nts = []
    for contig in relevant_contigs:
        block_str = contig.split('_')[1]
        site_positions = np.where(np.array(list(block_str)) == motif_ind)[0] * block_size + block_center
        for pos in site_positions:
            reference_motif = reference[contig][(pos - 2):(pos + 3)]
            this_tPos_nts = get_nucleotides_aligned_to_target_pos(bam_in, contig, pos, predStr_features, reference_motif, enforce_ref_5mer)
            if len(this_tPos_nts) > 0:
                motif_nts.extend(this_tPos_nts)

    return motif_nts

def get_nucleotides_aligned_to_target_pos(alignment, contig, target_pos, dict_predStr_feature, ref_motif=None, enforce_ref_5mer=False):
    all_nts = []
    for pileupcolumn in alignment.pileup(contig, target_pos, target_pos+1, truncate=True):
        if pileupcolumn.pos == target_pos:
            valid_counts = 0
            for ind, pileupread in enumerate(pileupcolumn.pileups):
                flag = pileupread.alignment.flag
                if flag!=0:
                    continue
                query_name = pileupread.alignment.query_name
                query_position = pileupread.query_position
                if query_position is None:
                    continue
                query_motif = pileupread.alignment.query_sequence[(query_position-2):(query_position+3)]
                if enforce_ref_5mer and (query_motif != ref_motif):
                    continue
                if query_name in dict_predStr_feature.keys():
                    this_read_predStr, this_read_feature = dict_predStr_feature[query_name]
                    this_site_motif = this_read_predStr[(query_position-2):(query_position+3)]
                    if this_site_motif!=query_motif:
                        print('!!! Error: Site motif {} =/= query {}!!!'.format(this_site_motif, query_motif))
                        print('Flag {}'.format(flag))
                        continue
                    this_site_feature = this_read_feature[query_position]
                    # motif_features_qPos[query_name] = (this_site_motif, this_site_feature, query_position)
                    all_nts.append(nucleotide(query_name, query_position, this_site_motif, this_site_feature))
                    valid_counts += 1
    return all_nts

########################################################################################################################
### mRNA ###############################################################################################################
########################################################################################################################
def get_nucleotides_from_multiple_reads(model, device, config, ext_layer, in_aligned_reads, batch_size=256):
    all_chunks = []
    chunk_sizes = []
    for this_aligned_read in in_aligned_reads:
        this_chunk = segment(this_aligned_read.norm_signal, config.seqlen)
        all_chunks.append(this_chunk)
        chunk_sizes.append(this_chunk.shape[0])
    if len(all_chunks)==0:
        return {}
    all_chunks = np.vstack(all_chunks)
    cum_chunk_sizes = np.cumsum([0] + chunk_sizes)

    event = torch.unsqueeze(torch.FloatTensor(all_chunks), 1).to(device, non_blocking=True)
    outs = []
    features = []
    for start in np.arange(0, event.shape[0], batch_size):
        stop = min(start+batch_size, event.shape[0])
        outs.append(model.forward(event[start:stop]))
        features.append(activation[ext_layer].detach().cpu().numpy())
    outs = torch.cat(outs, 1)
    features = np.vstack(features)

    all_nts = []
    for i, aligned_read in enumerate(in_aligned_reads):
        this_out = outs[:, cum_chunk_sizes[i]:cum_chunk_sizes[i+1], :]
        pred_str = []
        str_features = []
        for j in range(this_out.shape[1]):
            this_pred_label, pred_locs = ctc_decoder(torch.softmax(this_out[:, j, :], dim=-1).cpu().detach().numpy(), alphabet=alphabet)
            pred_str.extend(this_pred_label)
            str_features.append(features[cum_chunk_sizes[i]:cum_chunk_sizes[i + 1]][j][:, pred_locs].T)

        pred_str = ''.join(pred_str)[::-1]
        str_features = np.vstack(str_features)[::-1]

        pred_motif = ''.join(pred_str[aligned_read.pos-2:aligned_read.pos+3])

        if pred_motif!=aligned_read.query_5mer:
            print('\n!!! Error: Predicted motif {} =/= aligned {} !!!\n'.format(pred_motif, aligned_read.query_5mer))
            continue

        str_features = np.vstack(str_features)
        nt_feature = str_features[aligned_read.pos]

        all_nts.append(aligned_read.create_nucleotide(in_pred_5mer=pred_motif, in_feature=nt_feature))

    return all_nts

def get_nucleotides_aligned_to_site(model, device, config, ext_layer, alignment, index_read_ids, contig, site, strand, thresh_coverage=0, max_num_reads=1000, enforce_ref_5mer=False, ref_5mer='NNNNN'):
    all_aligned_reads = []
    for pileupcolumn in alignment.pileup(contig, site, site + 1, truncate=True):
        if pileupcolumn.pos == site:
            coverage = pileupcolumn.get_num_aligned()
            if coverage>thresh_coverage:
                valid_counts = 0
                for pileupread in pileupcolumn.pileups:
                    flag = pileupread.alignment.flag
                    if not (
                            ((strand=='+') and (flag==0))
                            or ((strand=='-') and (flag==16))
                    ):
                        continue
                    query_name = pileupread.alignment.query_name
                    query_position = pileupread.query_position
                    if query_position is None:
                        continue
                    if flag==16:
                        query_position = pileupread.alignment.query_length - query_position - 1
                    query_sequence = pileupread.alignment.get_forward_sequence()
                    query_5mer = query_sequence[(query_position - 2):(query_position + 3)]
                    # print('Strand {}, Flag {}, Ref 5mer: {}, Query 5-mer: {}'.format(strand, flag, ref_5mer, query_5mer))
                    if enforce_ref_5mer and (query_5mer!=ref_5mer):
                        continue
                    if (query_name in index_read_ids.keys()):
                        valid_counts += 1
                        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
                        all_aligned_reads.append(aligned_read(
                            read_id=query_name,
                            pos=query_position,
                            query_5mer=query_5mer,
                            norm_signal=this_read_signal,
                            flag=flag
                        ))
                    if (max_num_reads>0) and (len(all_aligned_reads)>=max_num_reads):
                        break
    if len(all_aligned_reads)==0:
        return []
    return get_nucleotides_from_multiple_reads(model, device, config, ext_layer, all_aligned_reads)