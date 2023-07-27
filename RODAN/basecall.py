#!/usr/bin/env python
#
# RODAN
# v1.0
# (c) 2020,2021,2022 Don Neumann
#

import torch
import numpy as np
import os, sys, argparse, time, glob
from fast_ctc_decode import beam_search, viterbi_search
from ont_fast5_api.fast5_interface import get_fast5_file
from torch.multiprocessing import Queue, Process
from models import Objectview, Rodan
import ont
import h5py
from tqdm import tqdm

vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }
alphabet = "".join(["N"] + list(vocab.values()))
alphabet_to_num = {v: k for k, v in enumerate(list(alphabet))}

default_rnaarch = [
    [-1, 256, 0, 3, 1, 1, 0],
    [-1, 256, 1, 10, 1, 1, 1],
    [-1, 256, 1, 10, 10, 1, 1],
    [-1, 320, 1, 10, 1, 1, 1],
    [-1, 384, 1, 15, 1, 1, 1],
    [-1, 448, 1, 20, 1, 1, 1],
    [-1, 512, 1, 25, 1, 1, 1],
    [-1, 512, 1, 30, 1, 1, 1],
    [-1, 512, 1, 35, 1, 1, 1],
    [-1, 512, 1, 40, 1, 1, 1],
    [-1, 512, 1, 45, 1, 1, 1],
    [-1, 512, 1, 50, 1, 1, 1],
    [-1, 768, 1, 55, 1, 1, 1],
    [-1, 768, 1, 60, 1, 1, 1],
    [-1, 768, 1, 65, 1, 1, 1],
    [-1, 768, 1, 70, 1, 1, 1],
    [-1, 768, 1, 75, 1, 1, 1],
    [-1, 768, 1, 80, 1, 1, 1],
    [-1, 768, 1, 85, 1, 1, 1],
    [-1, 768, 1, 90, 1, 1, 1],
    [-1, 768, 1, 95, 1, 1, 1],
    [-1, 768, 1, 100, 1, 1, 1]
]

def segment(seg, s):
    seg = np.concatenate((seg, np.zeros((-len(seg)%s))))
    nrows=((seg.size-s)//s)+1
    n=seg.strides[0]
    return np.lib.stride_tricks.as_strided(seg, shape=(nrows,s), strides=(s*n, n))


def convert_statedict(state_dict):
    from collections import OrderedDict
    new_checkpoint = OrderedDict()
    for k, v in state_dict.items():
        name = k[7:] # remove module.
        new_checkpoint[name] = v
    return new_checkpoint

activation = {}
def load_model(modelfile, config = None, args = None):
    if modelfile is None:
        sys.stderr.write("No model file specified!")
        sys.exit(1)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if args.debug: print("Using device:", device)
    if args.arch is not None:
        model = Rodan(config=config, arch=args.arch).to(device)
    else:
        model = Rodan(config=config).to(device)
    if args.debug: print("Loading pretrained weights:", modelfile)
    state_dict = torch.load(modelfile, map_location=device)["state_dict"]
    if "state_dict" in state_dict:
        model.load_state_dict(convert_statedict(state_dict["state_dict"]))
    else:
        model.load_state_dict(torch.load(modelfile, map_location=device)["state_dict"])
    if args.debug: print(model)

    def get_activation(name):
        def hook(model, input, output):
            activation[name] = output.detach()
        return hook

    model.get_submodule(args.extraction_layer).register_forward_hook(get_activation(args.extraction_layer))

    model.eval()
    torch.set_grad_enabled(False)

    return model, device


def mp_files(queue, config, args):
    dir = args.fast5dir
    if args.list_filenames is None:
        list_files = list(glob.iglob(dir+"/**/*.fast5", recursive=True))
    else:
        filenames = []
        with open(args.list_filenames, 'r') as f_in:
            list_files = [fname.rstrip('\n') for fname in f_in.readlines()]
        # list_files = [os.path.join(dir, path) for path in filenames]
    if args.debug: print('Running basecaller on files:\n', '\n'.join(list_files), flush=True)
    chunkname = []
    chunks = None
    queuechunks = None
    chunkremainder = None
    for file in list_files:
        try:
            f5 = get_fast5_file(file, mode="r")
        except:
            if args.debug: print('Error opening {}. Skipping...'.format(file), flush=True)
            continue
        for read in f5.get_reads():
            while queue.qsize() >= 100:
                time.sleep(1)
            #outfile = os.path.splitext(os.path.basename(file))[0]
            try:
                signal = read.get_raw_data(scale=True)
                if args.debug: print("mp_files:", file)
            except:
                continue
            signal_start = 0
            signal_end = len(signal)
            med, mad = ont.med_mad(signal[signal_start:signal_end])
            signal = (signal[signal_start:signal_end] - med) / mad
            newchunks = segment(signal, config.seqlen)
            if chunks is not None:
                chunks = np.concatenate((chunks, newchunks), axis=0)
                queuechunks += [read.read_id] * newchunks.shape[0]
            else:
                chunks = newchunks
                queuechunks = [read.read_id] * newchunks.shape[0]
            if chunks.shape[0] >= args.batchsize:
                for i in range(0, chunks.shape[0]//args.batchsize, args.batchsize):
                    queue.put((queuechunks[:args.batchsize], chunks[:args.batchsize]))
                    chunks = chunks[args.batchsize:]
                    queuechunks = queuechunks[args.batchsize:]
        f5.close()
    if len(queuechunks) > 0:
        if args.debug: print("queuechunks:", len(queuechunks), chunks.shape[0])
        for i in range(0, int(np.ceil(chunks.shape[0]/args.batchsize)), args.batchsize):
            start = i * args.batchsize
            end = start + args.batchsize
            if end > chunks.shape[0]: end = chunks.shape[0]
            queue.put((queuechunks[start:end], chunks[start:end]))
            if args.debug: print("put last chunk", chunks[start:end].shape[0])
    queue.put(("end", None))


def get_base_probs_and_activations(in_event, in_model, in_device):
    base_probs = in_model.forward(in_event)
    layer_activation = activation[args.extraction_layer]
    return base_probs, layer_activation


def mp_gpu(inqueue, outqueue, config, args):
    model, device = load_model(args.model, config, args)
    shtensor = None
    actensor = None
    if device.type=='cpu':
        pin_memory = False
    else:
        pin_memory = True
    while True:
        time1 = time.perf_counter()
        read = inqueue.get()
        file = read[0]
        if type(file) == str: 
            outqueue.put(("end", None))
            break
        chunks = read[1]
        for i in range(0, chunks.shape[0], config.batchsize):
            end = i+config.batchsize
            if end > chunks.shape[0]: end = chunks.shape[0]
            event = torch.unsqueeze(torch.FloatTensor(chunks[i:end]), 1).to(device, non_blocking=True)

            if args.dump_features:
                out, activations = get_base_probs_and_activations(event, model, device)
            else:
                out = model.forward(event)

            if shtensor is None:
                shtensor = torch.empty((out.shape), pin_memory=pin_memory, dtype=out.dtype)
            if out.shape[1] != shtensor.shape[1]:
                shtensor = torch.empty((out.shape), pin_memory=pin_memory, dtype=out.dtype)
            logitspre = shtensor.copy_(out).numpy()
            if args.debug: print("mp_gpu:", logitspre.shape)

            if args.dump_features:
                if actensor is None:
                    actensor = torch.empty((activations.shape), pin_memory=pin_memory, dtype=activations.dtype)
                if activations.shape[0] != actensor.shape[0]:
                    actensor = torch.empty((activations.shape), pin_memory=pin_memory, dtype=activations.dtype)
                np_activations = actensor.copy_(activations).numpy()
                np_activations = np_activations.transpose((2, 0, 1))

                outqueue.put((file, logitspre, np_activations))
                del activations
                del np_activations
            else:
                outqueue.put((file, logitspre))

            del out
            del logitspre


def get_basecall_and_features(in_base_probs, layer_activation=None, dump_features=False):
    if args.decoder == 'beam':
        decoder = beam_search
    elif args.decoder == 'viterbi':
        decoder = viterbi_search

    pred_labels = []
    out_features = []
    for i in range(in_base_probs.shape[1]):
        exp = np.exp(in_base_probs[:, i, :])
        probs = exp / np.sum(exp, axis=-1, keepdims=True)
        this_pred_label, pred_locs = decoder(probs, alphabet=alphabet)
        pred_labels.append(this_pred_label)

        if dump_features:
            num_locs = layer_activation.shape[0]
            pred_locs_corrected = []
            for loc_ind in range(len(this_pred_label)):
                start_loc = pred_locs[loc_ind]
                if loc_ind < (len(this_pred_label) - 1):
                    end_loc = pred_locs[loc_ind + 1]
                else:
                    end_loc = probs.shape[0]
                pred_locs_corrected.append(
                    start_loc + np.argmax(probs[start_loc:end_loc, alphabet_to_num[this_pred_label[loc_ind]]]))

            if args.feature_width == 0:
                this_feature = layer_activation[pred_locs_corrected, i, :]
            else:
                this_feature = []
                for loc_shift in range(-args.feature_width, args.feature_width + 1):
                    shifted_locs = [max(min(x + loc_shift, num_locs - 1), 0) for x in pred_locs_corrected]
                    this_feature.append(layer_activation[shifted_locs, i, :])
                this_feature = np.hstack(this_feature)
            out_features.append(this_feature)

    pred_labels = ''.join(pred_labels)

    if args.reverse:
        pred_labels = pred_labels[::-1]

    if dump_features:
        out_features = np.vstack(out_features)
        if args.reverse:
            out_features = out_features[::-1]
        return pred_labels, out_features
    else:
        return pred_labels, None


def mp_write(queue, config, args):
    files = None
    chunks = None
    totprocessed = 0
    finish = False

    read_features = {}
    with open(os.path.join(args.outdir, f'rodan.fasta'), 'w') as h_basecall:
        while True:
            if queue.qsize() > 0:
                newchunk = queue.get()
                if type(newchunk[0]) == str:
                    if not len(files): break
                    finish = True
                else:
                    if chunks is not None:
                        if args.dump_features:
                            activations = np.concatenate((activations, newchunk[2]), axis=1)
                        chunks = np.concatenate((chunks, newchunk[1]), axis=1)
                        files = files + newchunk[0]
                    else:
                        if args.dump_features:
                            activations = newchunk[2]
                        chunks = newchunk[1]
                        files = newchunk[0]

                while files.count(files[0]) < len(files) or finish:
                    totlen = files.count(files[0])
                    callchunk = chunks[:, :totlen, :]
                    if args.dump_features:
                        actichunk = activations[:, :totlen, :]
                    else:
                        actichunk = None
                    seq, features = get_basecall_and_features(callchunk, actichunk, args.dump_features)

                    readid = os.path.splitext(os.path.basename(files[0]))[0]
                    h_basecall.write(">" + readid + "\n")
                    h_basecall.write(seq + "\n")

                    newchunks = chunks[:, totlen:, :]
                    chunks = newchunks
                    files = files[totlen:]

                    if args.dump_features:
                        read_features[readid] = features
                        newactivations = activations[:, totlen:, :]
                        activations = newactivations

                    totprocessed += 1
                    if totprocessed%500==0: print(f'{totprocessed} reads processed', flush=True)
                    if finish and not len(files): break
                if finish: break
        print(f'Total {totprocessed} reads')

    if args.dump_features:
        print('Now dumping features...')
        with h5py.File(os.path.join(args.outdir, f'features.h5'), 'w') as h_features:
            for id, feat in tqdm(read_features.items()):
                h_features.create_dataset(id, data=feat)


if __name__ == "__main__":
    tic = time.time()

    parser = argparse.ArgumentParser(description='Basecall fast5 files')
    parser.add_argument("--fast5dir", default=None, type=str)
    parser.add_argument("--outdir", type=str)
    parser.add_argument("--list_filenames", default=None)
    parser.add_argument("-a", "--arch", default=None, type=str, help="architecture settings")
    parser.add_argument("-m", "--model", default="rna.torch", type=str, help="default: rna.torch")
    parser.add_argument("-r", "--reverse", default=True, action="store_true", help="reverse for RNA (default: True)")
    parser.add_argument("-b", "--batchsize", default=200, type=int, help="default: 200")
    parser.add_argument("--decoder", default="viterbi", help="beam_search or viterbi")
    parser.add_argument("--extraction_layer", default="convlayers.conv21")
    parser.add_argument('--feature_width', type=int, default=0)
    parser.add_argument("-B", "--beamsize", default=5, type=int, help="CTC beam search size (default: 5)")
    parser.add_argument("-e", "--errors", default=False, action="store_true")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--dump_features", default=False, action="store_true")
    args = parser.parse_args()

    print('Running RODAN basecaller', flush=True)

    os.makedirs(args.outdir, exist_ok=True)

    torchdict = torch.load(args.model, map_location="cpu")
    origconfig = torchdict["config"]

    if args.debug:
        print(origconfig)
    origconfig["debug"] = args.debug
    config = Objectview(origconfig)
    config.batchsize = args.batchsize

    if args.arch != None:
        if args.debug: print("Loading architecture from:", args.arch)
        args.arch = eval(open(args.arch, "r").read())
    else:
        args.arch = default_rnaarch

    if args.debug: print("Using sequence len:", int(config.seqlen))
    
    torch.backends.cudnn.enabled = True
    torch.backends.cudnn.deterministic = True

    call_queue = Queue()
    write_queue = Queue()
    p1 = Process(target=mp_files, args=(call_queue, config, args,))
    p2 = Process(target=mp_gpu, args=(call_queue, write_queue, config, args,))
    p3 = Process(target=mp_write, args=(write_queue, config, args))
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()

    toc = time.time()
    print('Finished in {:.1f} mins'.format((toc-tic)/60), flush=True)
