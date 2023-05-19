import numpy as np

class nucleotide:
    def __init__(self, read_id='', pos=-1, pred_5mer='NNNNN', feature=[], mod_prob=-1):
        self.read_id = str(read_id)
        self.pos = int(pos)
        self.pred_5mer = str(pred_5mer)
        self.feature = np.array(feature)
        self.mod_prob = float(mod_prob)

class aligned_read:
    def __init__(self, read_id='', pos=-1, query_5mer='NNNNN', pred_5mer='NNNNN', norm_signal=[], flag=-1):
        self.read_id = str(read_id)
        self.pos = int(pos)
        self.query_5mer = str(query_5mer)
        self.pred_5mer = str(pred_5mer)
        self.norm_signal = np.array(norm_signal)
        self.flag = int(flag)

    def create_nucleotide(self, in_pred_5mer, in_feature):
        return nucleotide(
            read_id = self.read_id,
            pos = self.pos,
            pred_5mer = in_pred_5mer,
            feature = in_feature,
        )