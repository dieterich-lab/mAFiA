import os


def get_read_ids(in_file):
    read_ids = []
    with open(in_file, 'r') as f_in:
        for this_line in f_in.readlines():
            read_ids.append(this_line.rstrip('\n'))
    return read_ids


raw_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/mouse_heart'
res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

# ds = 'HFpEF'
# raw_conds = ['LW_ctrl_H68_RTA', 'LW_ctrl_H69_RTA']
# res_cond = 'ctrl_merged'
# raw_conds = ['LW_HFpEF_H54_RTA', 'LW_HFpEF_H55_RTA']
# res_cond = 'HFpEF_merged'

# ds = 'Diet'
# raw_conds = ['A1_WT_CD', 'A2_WT_CD']
# res_cond = 'WT_CD_merged'
# raw_conds = ['C1_WT_WD', 'C2_WT_WD']
# res_cond = 'WT_WD_merged'

ds = 'TAC'
# raw_conds = ['40-31', '40-27', '40-28', '40-32']
# res_cond = 'SHAM_merged'
raw_conds = ['40-29', '40-26', '40-33', '40-30']
res_cond = 'TAC_merged'

res_read_ids = get_read_ids(os.path.join(res_dir, ds, res_cond, f'read_ids.txt'))
raw_read_ids = {this_cond: get_read_ids(os.path.join(raw_dir, ds, f'read_ids_{this_cond}.txt'))
                for this_cond in raw_conds}

print(f'res size: {len(res_read_ids)}')
for this_cond in raw_conds:
    num_intersection = len(set(res_read_ids).intersection(set(raw_read_ids[this_cond])))
    print(f'res intersect with {this_cond}: {num_intersection} / {len(raw_read_ids[this_cond])}')