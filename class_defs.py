import numpy as np

class nucleotide:
    def __init__(self, read_id='', pos=-1, pred_5mer='NNNNN', feature=[], mod_prob=-1):
        self.read_id = str(read_id)
        self.pos = int(pos)
        self.pred_5mer = str(pred_5mer)
        self.feature = np.array(feature)
        self.mod_prob = float(mod_prob)