import os
import requests, sys
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def get_genomic_coord_from_cDNA(id, pos):
    server = "http://rest.ensembl.org"
    ext = "/map/cdna/{}/{}..{}?".format(id, pos, pos+1)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    # print(repr(decoded))

    return decoded['mappings'][0]

test_dataset = '75_WT_25_IVT'
m6Anet_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Christoph/m6anet/workflow_tx/inference/{}/data.site_proba.csv'.format(test_dataset)
glori_file = '/home/adrian/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv'
out_dir = '/home/adrian/img_out/m6Anet'
os.makedirs(out_dir, exist_ok=True)

df_m6Anet = pd.read_csv(m6Anet_file)
df_glori = pd.read_csv(glori_file)

df_m6Anet_filtered = df_m6Anet
# df_m6Anet_filtered = df_m6Anet[df_m6Anet['probability_modified']>0.9]

dict_chr = {
    str(i) : 'chr{}'.format(i) for i in range(1, 23)
}
dict_chr['MT'] = 'chrM'

collected_sites = []
for ind, row in tqdm(df_m6Anet_filtered.iterrows()):
    # print(ind)
    tid = row['transcript_id'].split('.')[0]
    tpos = row['transcript_position']
    mod_ratio = row['mod_ratio']

    mapping = get_genomic_coord_from_cDNA(tid, tpos)
    if mapping['seq_region_name'] not in dict_chr.keys():
        continue
    chr_name = dict_chr[mapping['seq_region_name']]

    sel_row = df_glori[(df_glori['Chr']==chr_name) * (df_glori['Sites']==mapping['start'])]
    if len(sel_row)>0:
        # print(sel_row)
        collected_sites.append((ind, sel_row['Ratio'].values[0], mod_ratio))

plt.figure(figsize=(8, 8))
plt.plot([v[1] for v in collected_sites], [v[2] for v in collected_sites], '.')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('GLORI', fontsize=15)
plt.ylabel('m6Anet', fontsize=15)
plt.title(test_dataset, fontsize=20)
plt.savefig()