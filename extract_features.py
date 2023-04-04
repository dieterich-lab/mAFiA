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
from time import time

# EXTRACTION_LAYER = 'convlayers.conv21.sqex'

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
    # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # print("Using device:", device)
    device, mem = get_freer_device()
    print("Using device {} with free memory {}MB".format(device, mem), flush=True)
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

def get_features_from_collection_of_signals(model, device, config, index_read_ids, layer, feature_width=0):
    print('Extracting features from {}, width {}...'.format(layer, feature_width))
    id_predStr_feature = {}
    for query_name in tqdm(index_read_ids.keys()):
        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
        # tic = time()
        this_read_features, this_read_predStr = get_features_from_signal(model, device, config, this_read_signal, layer, feature_width)
        # toc = time()
        # print('Time: {:.3f}'.format(toc-tic))
        id_predStr_feature[query_name] = (this_read_predStr, this_read_features)
    return id_predStr_feature

# def extract_features_from_signal(model, device, config, signal, pos, bam_motif):
#     chunks = segment(signal, config.seqlen)
#     event = torch.unsqueeze(torch.FloatTensor(chunks), 1).to(device, non_blocking=True)
#     tic = time()
#     out = model.forward(event)
#     toc = time()
#     print('v1 model runtime: {:.2f}'.format(toc-tic))
#
#     pred_labels = []
#     features = []
#     for i in range(out.shape[1]):
#         this_pred_label, pred_locs = ctc_decoder(torch.softmax(out[:, i, :], dim=-1).cpu().detach().numpy(), alphabet=alphabet)
#         pred_labels.append(this_pred_label)
#         this_feature = activation['conv21'].detach().cpu().numpy()[i, :, pred_locs]
#         features.append(this_feature)
#
#     pred_label = ''.join(pred_labels)[::-1]
#     features = np.vstack(features)[::-1]
#
#     pred_motif = pred_label[pos-2:pos+3]
#     # print('Predicted motif {}, aligned {}, matching {}'.format(pred_motif, bam_motif, pred_motif == bam_motif))
#     if pred_motif!=bam_motif:
#         print('\n!!!Error: Predicted motif =/= aligned!!!\n')
#         return None, pred_motif
#
#     return features[pos], pred_motif

# def extract_features_from_signal_v2(model, device, config, signal, pos, bam_motif):
#     chunks = segment(signal, config.seqlen)
#     event = torch.unsqueeze(torch.FloatTensor(chunks), 1).to(device, non_blocking=True)
#     tic = time()
#     out = model.forward(event)
#     toc = time()
#     print('v2 model runtime: {:.2f}'.format(toc-tic))
#
#     pred_labels = []
#     features = None
#     for i in range(out.shape[1])[::-1]:
#         this_pred_label, pred_locs = ctc_decoder(torch.softmax(out[:, i, :], dim=-1).cpu().detach().numpy(), alphabet=alphabet)
#         pred_labels.extend(this_pred_label[::-1])
#         this_feature = activation['conv21'].detach().cpu().numpy()[i, :, pred_locs]
#         if features is None:
#             features = this_feature[::-1]
#         else:
#             features = np.vstack((features, this_feature[::-1]))
#         if len(pred_labels)>=(pos+3):
#             break
#
#     pred_label = ''.join(pred_labels)
#     # all_features = np.vstack(features)[::-1]
#
#     pred_motif = pred_label[pos-2:pos+3]
#     # print('Predicted motif {}, aligned {}, matching {}'.format(pred_motif, bam_motif, pred_motif == bam_motif))
#     if pred_motif!=bam_motif:
#         print('\n!!!!!!!!!!!!!!!!!!!Error: Predicted motif =/= aligned!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
#         return None, pred_motif
#
#     return features[pos], pred_motif

def extract_features_from_multiple_signals(model, device, config, ext_layer, site_normReads_qPos_motif_sense, batch_size=128):
    qNames = []
    all_chunks = []
    chunk_sizes = []
    all_pos = []
    all_motifs = []
    all_senses = []
    for qname, (this_signal, this_pos, this_motif, this_sense) in site_normReads_qPos_motif_sense.items():
        qNames.append(qname)
        this_chunk = segment(this_signal, config.seqlen)
        all_chunks.append(this_chunk)
        chunk_sizes.append(this_chunk.shape[0])
        all_pos.append(this_pos)
        all_motifs.append(this_motif)
        all_senses.append(this_sense)
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

    site_motif_features = {}
    for i in range(len(all_pos)):
        this_qname = qNames[i]
        this_pos = all_pos[i]
        this_sense = all_senses[i]
        this_out = outs[:, cum_chunk_sizes[i]:cum_chunk_sizes[i+1], :]
        pred_str = []
        str_features = []
        for j in range(this_out.shape[1]):
            this_pred_label, pred_locs = ctc_decoder(torch.softmax(this_out[:, j, :], dim=-1).cpu().detach().numpy(), alphabet=alphabet)
            pred_str.extend(this_pred_label)
            str_features.append(features[cum_chunk_sizes[i]:cum_chunk_sizes[i + 1]][j][:, pred_locs].T)

        pred_str = ''.join(pred_str)[::-1]
        str_features = np.vstack(str_features)[::-1]

        query_motif = all_motifs[i]
        if this_sense=='+':
            pred_motif = ''.join(pred_str[this_pos-2:this_pos+3])
        else:
            pred_motif = ''.join(pred_str[::-1][this_pos-2:this_pos+3][::-1])

        if pred_motif!=query_motif:
            print('\n!!! Error: Predicted motif {} =/= aligned {} !!!\n'.format(pred_motif, query_motif))
            continue

        str_features = np.vstack(str_features)
        if this_sense=='+':
            site_motif_features[this_qname] = (pred_motif, str_features[this_pos])
        else:
            site_motif_features[this_qname] = (pred_motif, str_features[::-1][this_pos])

    return site_motif_features


def collect_all_motif_features(motif_ind, reference, bam_in, predStr_features, enforce_motif=False):
    BLOCK_CENTER = 16
    BLOCK_SIZE = 33

    relevant_contigs = [k for k in reference.keys() if motif_ind in k.split('_')[1]]
    motif_features = []
    for contig in relevant_contigs:
        block_str = contig.split('_')[1]
        site_positions = np.where(np.array(list(block_str)) == motif_ind)[0] * BLOCK_SIZE + BLOCK_CENTER
        for pos in site_positions:
            reference_motif = reference[contig][(pos - 2):(pos + 3)]
            if enforce_motif:
                site_motif_features = collect_site_features(bam_in, contig, pos, predStr_features, reference_motif)
            else:
                site_motif_features = collect_site_features(bam_in, contig, pos, predStr_features)
            if len(site_motif_features) > 0:
                motif_features.append(site_motif_features)
    return {k: v for dic in motif_features for k, v in dic.items()}

def collect_site_features(alignment, contig, pos, dict_predStr_feature, enforce_motif=None):
    site_motif_features = {}
    for pileupcolumn in alignment.pileup(contig, pos, pos+1, truncate=True):
        if pileupcolumn.pos == pos:
            valid_counts = 0
            for ind, pileupread in enumerate(pileupcolumn.pileups):
                query_name = pileupread.alignment.query_name
                query_position = pileupread.query_position
                if query_position is None:
                    continue
                query_motif = pileupread.alignment.query_sequence[(query_position-2):(query_position+3)]
                flag = pileupread.alignment.flag
                if (enforce_motif is not None) and (query_motif != enforce_motif):
                    continue
                if (query_position and (flag == 0) and (query_name in dict_predStr_feature.keys())):
                    this_read_predStr, this_read_feature = dict_predStr_feature[query_name]
                    this_site_motif = this_read_predStr[(query_position-2):(query_position+3)]
                    if this_site_motif!=query_motif:
                        print('!!! Error: Site motif {} =/= query {}!!!'.format(this_site_motif, query_motif))
                        continue
                    this_site_feature = this_read_feature[query_position]
                    site_motif_features[query_name] = (this_site_motif, this_site_feature)
                    valid_counts += 1
    return site_motif_features

# def collect_features_from_aligned_site(model, device, config, alignment, index_read_ids, contig, site, thresh_coverage=0, enforce_motif=None):
#     MAX_READS_IN_PILEUP = 500
#     site_motif_features = {}
#     for pileupcolumn in alignment.pileup(contig, site, site + 1, truncate=True, min_base_quality=20, min_mapping_quality=20):
#         if pileupcolumn.pos == site:
#             coverage = pileupcolumn.get_num_aligned()
#             print('Coverage {}'.format(coverage), flush=True)
#             if coverage>thresh_coverage:
#                 valid_counts = 0
#                 for ind, pileupread in enumerate(pileupcolumn.pileups):
#                     if valid_counts >= MAX_READS_IN_PILEUP:
#                         break
#                     # print('Pileup {}/{}'.format(ind, coverage))
#                     query_name = pileupread.alignment.query_name
#                     # query_position = pileupread.query_position_or_next
#                     query_position = pileupread.query_position
#                     if query_position is None:
#                         continue
#                     query_motif = pileupread.alignment.query_sequence[query_position - 2:query_position + 3]
#                     flag = pileupread.alignment.flag
#                     if (enforce_motif is not None) and (query_motif!=enforce_motif):
#                         continue
#                     if query_position and (flag==0) and (query_name in index_read_ids.keys()):
#                         valid_counts += 1
#                         query_motif = pileupread.alignment.query_sequence[query_position-2:query_position+3]
#                         this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
#                         # this_read_signal = id_signal[query_name]
#                         tic = time()
#                         this_read_features_v1, this_read_motif_v1 = extract_features_from_signal(model, device, config, this_read_signal, query_position, query_motif)
#                         toc = time()
#                         print('v1 time: {:.3f}'.format(toc-tic))
#                         tic = time()
#                         this_read_features, this_read_motif = extract_features_from_signal_v2(model, device, config, this_read_signal, query_position, query_motif)
#                         toc = time()
#                         print('v2 time: {:.3f}'.format(toc-tic))
#                         print(this_read_motif_v1, this_read_motif)
#                         print(np.all(this_read_features_v1==this_read_features))
#                         if this_read_features is not None:
#                             site_motif_features[query_name] = (this_read_motif, this_read_features)
#     return site_motif_features

def collect_features_from_aligned_site_v2(model, device, config, ext_layer, alignment, index_read_ids, contig, site, thresh_coverage=0, max_num_reads=1000, enforce_motif=None):
    site_normReads_qPos_motif_sense = {}
    for pileupcolumn in alignment.pileup(contig, site, site + 1, truncate=True):
        if pileupcolumn.pos == site:
            coverage = pileupcolumn.get_num_aligned()
            if coverage>thresh_coverage:
                valid_counts = 0
                for pileupread in pileupcolumn.pileups:
                    query_name = pileupread.alignment.query_name
                    # query_position = pileupread.query_position_or_next
                    query_position = pileupread.query_position
                    if query_position is None:
                        continue
                    query_motif = pileupread.alignment.query_sequence[(query_position-2):(query_position+3)]
                    flag = pileupread.alignment.flag
                    if flag==16:
                        query_motif = str(Seq(query_motif).reverse_complement())
                    if (enforce_motif is not None) and (query_motif!=enforce_motif):
                        continue
                    if query_position and (flag in [0, 16]) and (query_name in index_read_ids.keys()):
                        valid_counts += 1
                        this_read_signal = get_norm_signal_from_read_id(query_name, index_read_ids)
                        # this_read_signal = id_signal[query_name]
                        sense = '+' if flag==0 else '-'
                        site_normReads_qPos_motif_sense[query_name] = (this_read_signal, query_position, query_motif, sense)
                    if len(site_normReads_qPos_motif_sense) >= max_num_reads:
                        break
    if len(site_normReads_qPos_motif_sense)==0:
        return {}
    site_motif_features = extract_features_from_multiple_signals(model, device, config, ext_layer, site_normReads_qPos_motif_sense)
    return site_motif_features