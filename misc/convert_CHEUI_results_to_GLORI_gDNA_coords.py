import os
from os.path import expanduser
HOME = expanduser('~')
import requests, sys
import pandas as pd
import numpy as np
from tqdm import tqdm

test_dataset = '100_WT_0_IVT'

# PRJ = os.path.join(HOME, 'Data')
# PRJ = '/prj'

def get_genomic_coord_from_cDNA(id, pos):
    server = "http://rest.ensembl.org"
    ext = "/map/cdna/{}/{}..{}?".format(id, pos, pos+1)
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    except:
        # print('Request failed')
        return None
    if not r.ok:
        # r.raise_for_status()
        # sys.exit()
        # print('Request failed')
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

result_dir = os.path.join('/scratch/achan/CHEUI')
cheui_file = os.path.join(result_dir, f'{test_dataset}/site_level_m6A_predictions.txt')
tx_file = os.path.join(HOME, 'Data/transcriptomes/GRCh38_102.bed')
glori_file = os.path.join(HOME, 'Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
out_dir = result_dir
os.makedirs(out_dir, exist_ok=True)

df_cheui = pd.read_csv(cheui_file, sep='\t')
df_tx = pd.read_csv(tx_file, sep='\t', usecols=list(range(4))+[5]+list(range(9, 12)),
                    names=['chr', 'start', 'stop', 'transcript_id', 'strand', 'num_blocks', 'block_sizes', 'block_starts'],
                    dtype={'chr': str})
df_glori = pd.read_csv(glori_file)

df_cheui_filtered = df_cheui

dict_chr = {
    str(i) : 'chr{}'.format(i) for i in range(1, 23)
}
dict_chr['MT'] = 'chrM'
dict_chr['X'] = 'chrX'
dict_chr['Y'] = 'chrY'

collected_sites = []
for ind, row in tqdm(df_cheui_filtered.iterrows()):
    # print(ind)
    tid = row['contig'].split('.')[0]
    tpos = row['position']
    mod_ratio = row['stoichiometry']

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
print('{} GLORI sites collected'.format(len(collected_sites)))

### slice df and check against ensembl API ###
df_cheui_glori = df_cheui.loc[[site[0] for site in collected_sites]]
df_cheui_glori['Chr'] = [site[1] for site in collected_sites]
df_cheui_glori['Sites'] = [site[2] for site in collected_sites]
df_cheui_glori['Ratio'] = [site[3] for site in collected_sites]
df_cheui_glori['Pvalue'] = [site[5] for site in collected_sites]
df_cheui_glori.to_csv(cheui_file.replace('.txt', '_glori.txt'), sep='\t')

bad_indices = []
for ind, row in tqdm(df_cheui_glori.iterrows()):
    # print(ind)
    ### ensembl API ###
    tid = row['contig'].split('.')[0]
    tpos = row['position']
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
        # print('{}: ensembl gpos {} != my gpos {}\n'.format(ind, ens_start, row['Sites']))
    # else:
    #     print('Okay\n')
print('{} bad cDNA->gDNA conversions out of {}'.format(len(bad_indices), len(df_cheui_glori)))

df_cheui_glori_filtered = df_cheui_glori.drop(bad_indices)
df_cheui_glori_filtered.to_csv(cheui_file.replace('.txt', '.glori_filtered.txt'), sep='\t')