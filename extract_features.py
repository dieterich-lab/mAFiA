import os
HOME = os.path.expanduser('~')
from utils import segment
import torch
from models import objectview

model_path = os.path.join(HOME, 'pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch')
torchdict = torch.load(model_path, map_location="cpu")
origconfig = torchdict["config"]
config = objectview(origconfig)
config.batchsize = 16

def extract_features_from_signal(signal, pos):
    chunks = segment(signal, config.seqlen)
    return