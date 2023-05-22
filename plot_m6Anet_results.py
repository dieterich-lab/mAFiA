import os
from os.path import expanduser
HOME = expanduser('~')
import requests, sys
import pandas as pd
import numpy as np
from tqdm import tqdm
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

PRJ = os.path.join(HOME, 'Data')
# PRJ = '/prj'

def get_genomic_coord_from_cDNA(id, pos):
    server = "http://rest.ensembl.org"
    ext = "/map/cdna/{}/{}..{}?".format(id, pos, pos+1)
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    except:
        print('Request failed')
        return None
    if not r.ok:
        # r.raise_for_status()
        # sys.exit()
        print('Request failed')
        return None
    decoded = r.json()
    # print(repr(decoded))

    return decoded['mappings'][0]

def get_cDNA_coords(in_tid, in_tpos, tx):
    sub_tx = tx[tx['transcript_id'] == in_tid]
    if len(sub_tx)==0:
        return None, None
    out_chr = sub_tx['chr'].values[0]
    num_blocks = sub_tx['num_blocks'].values[0]
    block_sizes = np.array([int(size) for size in sub_tx['block_sizes'].values[0].split(',') if len(size)>0])
    block_starts = np.array([int(size) for size in sub_tx['block_starts'].values[0].split(',') if len(size)>0])

    strand = sub_tx['strand'].values[0]
    # print('{}, strand {}, {} blocks'.format(in_tid, strand, num_blocks))
    if strand=='+':
        block_cumsum = np.concatenate([[0], np.cumsum(block_sizes)])
        block_ind = np.where(in_tpos >= block_cumsum)[0][-1]
        local_gpos = block_starts[block_ind] + in_tpos - block_cumsum[block_ind]
    else:
        block_cumsum = np.concatenate([[0], np.cumsum(block_sizes[::-1])])
        block_ind = np.where(in_tpos >= block_cumsum)[0][-1]
        local_gpos = block_starts[::-1][block_ind] + block_sizes[::-1][block_ind] - (in_tpos - block_cumsum[block_ind])

    out_gpos = sub_tx['start'].values[0] + local_gpos

    return out_chr, out_gpos

test_dataset = '0_WT_100_IVT'
m6Anet_file = os.path.join(PRJ, 'TRR319_RMaP/Project_BaseCalling/Christoph/m6anet/workflow_tx/inference/{}/data.site_proba.csv'.format(test_dataset))
tx_file = os.path.join(HOME, 'Data/transcriptomes/GRCh38_102.bed')
glori_file = os.path.join(HOME, 'Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
out_dir = os.path.join(HOME, 'img_out/m6Anet')
os.makedirs(out_dir, exist_ok=True)

df_m6Anet = pd.read_csv(m6Anet_file)
df_tx = pd.read_csv(tx_file, sep='\t', usecols=list(range(4))+[5]+list(range(9, 12)), names=['chr', 'start', 'stop', 'transcript_id', 'strand', 'num_blocks', 'block_sizes', 'block_starts'])
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

    ### parsing bed file ###
    contig, gpos = get_cDNA_coords(tid, tpos, df_tx)
    if ((contig is None) or (contig not in dict_chr.keys())):
        continue

    ### check ###
    # print('contig {}, pos {}'.format(contig==mapping['seq_region_name'], gpos==mapping['start']))

    sel_row = df_glori[(df_glori['Chr']==dict_chr[contig]) * (df_glori['Sites']==gpos)]
    if len(sel_row)>0:
        # print(sel_row)
        collected_sites.append((ind, dict_chr[contig], gpos, sel_row['Ratio'].values[0], mod_ratio, sel_row['Pvalue'].values[0]))

### slice df and check against ensembl API ###
df_m6Anet_glori = df_m6Anet.loc[[site[0] for site in collected_sites]]
df_m6Anet_glori['Chr'] = [site[1] for site in collected_sites]
df_m6Anet_glori['Sites'] = [site[2] for site in collected_sites]
df_m6Anet_glori['GLORI'] = [site[3] for site in collected_sites]
df_m6Anet_glori.to_csv(m6Anet_file.replace('.csv', '_glori.csv'))

bad_indices = []
for ind, row in tqdm(df_m6Anet_glori.iterrows()):
    # print(ind)
    ### ensembl API ###
    tid = row['transcript_id'].split('.')[0]
    tpos = row['transcript_position']
    mapping = get_genomic_coord_from_cDNA(tid, tpos)
    if mapping is None:
        bad_indices.append(ind)
        continue

    if mapping['seq_region_name'] not in dict_chr.keys():
        continue
    ens_chr = dict_chr[mapping['seq_region_name']]
    ens_start = mapping['start']
    if ens_start!=row['Sites']:
        bad_indices.append(ind)
        print('{}: ensembl gpos {} != my gpos {}\n'.format(ind, ens_start, row['Sites']))
    # else:
    #     print('Okay\n')

df_m6Anet_glori_filtered = df_m6Anet_glori.drop(bad_indices)
df_m6Anet_glori_filtered.to_csv(m6Anet_file.replace('.csv', '_glori_filtered.csv'))

### plot ###
thresh_n_reads = 50
thresh_pval = 1E-99
thresh_df = df_m6Anet_glori_filtered[
    (df_m6Anet_glori_filtered['n_reads']>=thresh_n_reads)
    * (df_m6Anet_glori_filtered['Pvalue']<thresh_pval)
    ]
plt.figure(figsize=(8, 8))
plt.plot(thresh_df['GLORI'], thresh_df['mod_ratio'], '.')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('GLORI', fontsize=15)
plt.ylabel('m6Anet', fontsize=15)
plt.title('{}\nn_reads$\\geq${}, p_val$\\leq${}'.format(test_dataset, thresh_n_reads, thresh_pval), fontsize=20)
plt.savefig(os.path.join(out_dir, '{}_nreads{}_pval{}.png'.format(test_dataset, thresh_n_reads, thresh_pval)), bbox_inches='tight')