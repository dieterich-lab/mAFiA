import os, sys
HOME = os.path.expanduser('~')
import torch
import numpy as np
from models import Objectview, Rodan
from fast_ctc_decode import beam_search, viterbi_search

CTC_MODE='viterbi'
if CTC_MODE=='beam':
    ctc_decoder = beam_search
elif CTC_MODE=='viterbi':
    ctc_decoder = viterbi_search

vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }
alphabet = "".join(["N"] + list(vocab.values()))
alphabet_to_num = {v: k for k, v in enumerate(list(alphabet))}

def convert_statedict(state_dict):
    from collections import OrderedDict
    new_checkpoint = OrderedDict()
    for k, v in state_dict.items():
        name = k[7:] # remove module.
        new_checkpoint[name] = v
    return new_checkpoint


class Backbone_Network:
    def __init__(self, model_path, extraction_layer, feature_width, batchsize=1024):
        print('Finding my backbone...')
        torchdict = torch.load(model_path, map_location="cpu")
        self.config = Objectview(torchdict["config"])
        self.extraction_layer = extraction_layer
        self.feature_width = feature_width
        self.batchsize = batchsize
        self.activation = {}
        self._load_model(model_path)
        print(f'Using device {self.device}, model {os.path.basename(model_path)} at extraction layer {extraction_layer}')

    def _segment(self, seg, s):
        seg = np.concatenate((seg, np.zeros((-len(seg) % s))))
        nrows = ((seg.size - s) // s) + 1
        n = seg.strides[0]
        return np.lib.stride_tricks.as_strided(seg, shape=(nrows, s), strides=(s * n, n))

    def _load_model(self, modelfile):
        if modelfile == None:
            sys.stderr.write("No model file specified!")
            sys.exit(1)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = Rodan(config=self.config).to(device)
        state_dict = torch.load(modelfile, map_location=device)["state_dict"]
        if "state_dict" in state_dict:
            model.load_state_dict(convert_statedict(state_dict["state_dict"]))
        else:
            model.load_state_dict(torch.load(modelfile, map_location=device)["state_dict"])
        model.get_submodule(self.extraction_layer).register_forward_hook(self._get_activation(self.extraction_layer))
        model.eval()
        torch.set_grad_enabled(False)

        self.model, self.device = model, device

    def _get_activation(self, name):
        def hook(model, input, output):
            self.activation[name] = output.detach()
        return hook

    def _get_base_probs_and_activations(self, in_chunks):
        event = torch.unsqueeze(torch.FloatTensor(in_chunks), 1).to(self.device, non_blocking=True)
        event_size = event.shape[0]
        if event_size <= self.batchsize:
            out = self.model.forward(event)
            layer_activation = self.activation[self.extraction_layer].detach().cpu().numpy()
        else:
            break_pts = np.arange(0, event_size, self.batchsize)
            start_stop_pts = [(start, stop) for start, stop in zip(break_pts, list(break_pts[1:]) + [event_size])]
            batch_out = []
            batch_layer_activation = []
            for (start, stop) in start_stop_pts:
                batch_out.append(self.model.forward(event[start:stop]))
                batch_layer_activation.append(self.activation[self.extraction_layer].detach().cpu().numpy())
            out = torch.concat(batch_out, 1)
            layer_activation = np.concatenate(batch_layer_activation, axis=0)
        return out, layer_activation


    def _get_basecall_and_features(self, in_base_probs, layer_activation):
        num_locs = layer_activation.shape[-1]
        pred_labels = []
        out_features = []
        for i in range(in_base_probs.shape[1]):
            probs = torch.softmax(in_base_probs[:, i, :], dim=-1).cpu().detach().numpy()
            this_pred_label, pred_locs = ctc_decoder(probs, alphabet=alphabet)

            pred_labels.append(this_pred_label)

            pred_locs_corrected = []
            for loc_ind in range(len(this_pred_label)):
                start_loc = pred_locs[loc_ind]
                if loc_ind < (len(this_pred_label) - 1):
                    end_loc = pred_locs[loc_ind + 1]
                else:
                    end_loc = probs.shape[0]
                pred_locs_corrected.append(
                    start_loc + np.argmax(probs[start_loc:end_loc, alphabet_to_num[this_pred_label[loc_ind]]]))

            if self.feature_width == 0:
                this_feature = layer_activation[i, :, pred_locs_corrected]
            else:
                this_feature = []
                for loc_shift in range(-self.feature_width, self.feature_width + 1):
                    shifted_locs = [max(min(x + loc_shift, num_locs - 1), 0) for x in pred_locs_corrected]
                    this_feature.append(layer_activation[i, :, shifted_locs])
                this_feature = np.hstack(this_feature)
            out_features.append(this_feature)
        pred_label = ''.join(pred_labels)[::-1]
        out_features = np.vstack(out_features)[::-1]

        return pred_label, out_features

    def get_features_from_signal(self, signal):
        chunks = self._segment(signal, self.config.seqlen)

        base_probs, activations = self._get_base_probs_and_activations(chunks)
        basecalls, features = self._get_basecall_and_features(base_probs, activations)

        return features, basecalls

    def _get_chunks_and_sizes_from_multiple_aligned_reads(self, in_aligned_reads):
        out_chunks = []
        chunk_sizes = []
        for this_aligned_read in in_aligned_reads:
            this_chunk = self._segment(this_aligned_read.norm_signal, self.config.seqlen)
            out_chunks.append(this_chunk)
            chunk_sizes.append(this_chunk.shape[0])
        if len(out_chunks) == 0:
            return [], -1

        out_chunks = np.vstack(out_chunks)
        cum_chunk_sizes = np.cumsum([0] + chunk_sizes)

        return out_chunks, cum_chunk_sizes

    def get_nucleotides_from_multiple_reads(self, in_aligned_reads):
        chunks, chunk_sizes = self._get_chunks_and_sizes_from_multiple_aligned_reads(in_aligned_reads)
        if len(chunks)==0:
            return []

        base_probs, activations = self._get_base_probs_and_activations(chunks)

        all_nts = []
        for i, aligned_read in enumerate(in_aligned_reads):
            this_chunk_base_probs = base_probs[:, chunk_sizes[i]:chunk_sizes[i + 1], :]
            this_chunk_activations = activations[chunk_sizes[i]:chunk_sizes[i + 1], :, :]
            this_chunk_basecalls, this_chunk_features = self._get_basecall_and_features(this_chunk_base_probs, this_chunk_activations)

            pred_motif = this_chunk_basecalls[aligned_read.read_pos-2 : aligned_read.read_pos+3]
            if pred_motif != aligned_read.query_5mer:
                print(f'\n!!! Error: Predicted motif {pred_motif} =/= aligned {aligned_read.query_5mer} !!!\n')
                continue

            nt_feature = this_chunk_features[aligned_read.read_pos]
            all_nts.append(aligned_read.create_nucleotide(in_pred_5mer=pred_motif, in_feature=nt_feature))

        return all_nts