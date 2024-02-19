import os
from os.path import expanduser
HOME = expanduser('~')
import requests, sys
import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq

# test_dataset = '100_WT_0_IVT'

# PRJ = os.path.join(HOME, 'Data')
# PRJ = '/prj'

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, 'fasta'):
    # if (record.id.isnumeric()) or (record.id in ['X', 'Y', 'MT']):
    ref[record.id] = record.seq

# cheui_file = sys.argv[1]
cheui_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/CHEUI/HEK293/100_WT_0_IVT/site_level_m6A_predictions.txt'

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

# def get_cDNA_coords(in_tid, in_tpos, tx):
#     sub_tx = tx[tx['transcript_id'] == in_tid]
#     if len(sub_tx)==0:
#         return None, None
#     out_chr = sub_tx['chr'].values[0]
#     num_blocks = sub_tx['num_blocks'].values[0]
#     block_sizes = np.array([int(size) for size in sub_tx['block_sizes'].values[0].split(',') if len(size)>0])
#     block_starts = np.array([int(size) for size in sub_tx['block_starts'].values[0].split(',') if len(size)>0])
#
#     strand = sub_tx['strand'].values[0]
#     # print('{}, strand {}, {} blocks'.format(in_tid, strand, num_blocks))
#     if strand=='+':
#         block_cumsum = np.concatenate([[0], np.cumsum(block_sizes)])
#         block_ind = np.where(in_tpos >= block_cumsum)[0][-1]
#         local_gpos = block_starts[block_ind] + in_tpos - block_cumsum[block_ind]
#     else:
#         block_cumsum = np.concatenate([[0], np.cumsum(block_sizes[::-1])])
#         block_ind = np.where(in_tpos >= block_cumsum)[0][-1]
#         local_gpos = block_starts[::-1][block_ind] + block_sizes[::-1][block_ind] - (in_tpos - block_cumsum[block_ind])
#
#     out_gpos = sub_tx['start'].values[0] + local_gpos
#
#     return out_chr, out_gpos

def combine_transcript_sites(df_in):
    chrom_pos = [tuple(pair) for pair in df_in[['chrom', 'chromStart']].values]
    collected_sites = []
    for this_chrom_pos in tqdm(set(chrom_pos)):
        sub_df = df_in[(df_in['chrom']==this_chrom_pos[0]) * df_in['chromStart']==this_chrom_pos[1]]
        strand = sub_df['strand'].values[0]
        if len(sub_df)>1:
            n_reads = sub_df['coverage'].values
            mod_ratio = sub_df['stoichiometry'].values
            coverage = np.sum(n_reads)
            avg_modRatio = np.sum(n_reads * mod_ratio) / coverage
        else:
            coverage = sub_df['coverage'].values[0]
            avg_modRatio = sub_df['stoichiometry'].values[0]

        ref5mer = sub_df['ref5mer'].values[0]
        pos = this_chrom_pos[1]
        if strand=='-':
            pos -= 1

        check5mer = ref[this_chrom_pos[0]][pos-2:pos+3]
        if strand=='-':
            check5mer = check5mer.reverse_complement()
        # print(strand, ref5mer, check5mer)
        if ref5mer == check5mer:
            collected_sites.append((this_chrom_pos[0], pos, strand, coverage, avg_modRatio, ref5mer))
    df_collected_sites = pd.DataFrame(collected_sites, columns=['chrom', 'chromStart', 'strand', 'coverage', 'modRatio', 'ref5mer'])
    df_collected_sites['chromEnd'] = df_collected_sites['chromStart'] + 1
    df_collected_sites = df_collected_sites[['chrom', 'chromStart', 'chromEnd', 'strand', 'coverage', 'modRatio', 'ref5mer']]
    df_collected_sites = df_collected_sites.sort_values(by=['chrom', 'chromStart'], ignore_index=True)
    df_collected_sites['modRatio'] = np.int32(np.round(df_collected_sites['modRatio'] * 100.0, 0))
    return df_collected_sites


# result_dir = os.path.join('/scratch/achan/CHEUI')
# cheui_file = os.path.join(result_dir, f'{test_dataset}/site_level_m6A_predictions.txt')
# tx_file = os.path.join(HOME, 'Data/transcriptomes/GRCh38_102.bed')
# glori_file = os.path.join(HOME, 'Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv')
# out_dir = result_dir
# os.makedirs(out_dir, exist_ok=True)

df_cheui = pd.read_csv(cheui_file, sep='\t')
# df_tx = pd.read_csv(tx_file, sep='\t', usecols=list(range(4))+[5]+list(range(9, 12)),
#                     names=['chr', 'start', 'stop', 'transcript_id', 'strand', 'num_blocks', 'block_sizes', 'block_starts'],
#                     dtype={'chr': str})
# df_glori = pd.read_csv(glori_file)

sel_motifs = [
    'GGACT',
    'GAACT',
    'GGACA',
    'AGACT',
    'GGACC',
    'TGACT'
]
df_cheui['ref5mer'] = [this_site[2:7]for this_site in df_cheui['site']]
df_cheui_filtered = df_cheui[df_cheui['ref5mer'].isin(sel_motifs)]

# df_cheui_filtered['position'] = df_cheui_filtered['position'] + 4

dict_chr = {
    str(i) : 'chr{}'.format(i) for i in range(1, 23)
}
# dict_chr['MT'] = 'chrM'
dict_chr['X'] = 'chrX'
# dict_chr['Y'] = 'chrY'

collected_sites = []
for ind, row in tqdm(df_cheui_filtered.iterrows()):
    # print(ind)
    tid = row['contig'].split('.')[0]
    tpos = row['position'] + 4
    mod_ratio = row['stoichiometry']

    ### parsing bed file ###
    # contig, gpos = get_cDNA_coords(tid, tpos, df_tx)
    # if ((contig is None) or (contig not in dict_chr.keys())):
    #     continue

    mapping = get_genomic_coord_from_cDNA(tid, tpos)
    if mapping is None:
        continue
    if mapping['seq_region_name'] not in dict_chr.keys():
        continue
    ens_chr = dict_chr[mapping['seq_region_name']]
    ens_start = mapping['start']

    # row['chrom'] = contig
    # row['chromStart'] = gpos

    row['chrom'] = ens_chr.lstrip('chr')
    row['chromStart'] = ens_start
    if mapping['strand'] == 1:
        row['strand'] = '+'
    else:
        row['strand'] = '-'

    collected_sites.append(row)

    ### check ###
    # print('contig {}, pos {}'.format(contig==mapping['seq_region_name'], gpos==mapping['start']))

    # sel_row = df_glori[(df_glori['Chr']==dict_chr[contig]) & (df_glori['Sites']==gpos)]
    # if len(sel_row)>0:
    #     # print(sel_row)
    #     collected_sites.append((ind, dict_chr[contig], gpos, sel_row['Ratio'].values[0], mod_ratio, sel_row['P_adjust'].values[0], sel_row['Strand'].values[0]))
print('{} sites collected'.format(len(collected_sites)))

df_out = pd.DataFrame(collected_sites)
df_out.to_csv(cheui_file+'.genome.6motifs', sep='\t', index=False)

df_out_combined = combine_transcript_sites(df_out)
df_out_combined.to_csv(cheui_file+'.genome.6motifs.combined', sep='\t', index=False)

### slice df and check against ensembl API ###
# df_cheui_glori = df_cheui.loc[[site[0] for site in collected_sites]]
# df_cheui_glori['Chr'] = [site[1] for site in collected_sites]
# df_cheui_glori['Sites'] = [site[2] for site in collected_sites]
# df_cheui_glori['Ratio'] = [site[3] for site in collected_sites]
# df_cheui_glori['P_adjust'] = [site[5] for site in collected_sites]
# df_cheui_glori['Strand'] = [str(site[6]) for site in collected_sites]
# df_cheui_glori.to_csv(cheui_file+'.glori', sep='\t', index=False)

# bad_indices = []
# for ind, row in tqdm(df_cheui_glori.iterrows()):
#     # print(ind)
#     ### ensembl API ###
#     tid = row['contig'].split('.')[0]
#     tpos = row['position'] + 4
#     mapping = get_genomic_coord_from_cDNA(tid, tpos)
#     if mapping is None:
#         bad_indices.append(ind)
#         continue
#
#     if mapping['seq_region_name'] not in dict_chr.keys():
#         continue
#     ens_chr = dict_chr[mapping['seq_region_name']]
#     ens_start = mapping['start']
#     if ens_start!=row['Sites']:
#         bad_indices.append(ind)
#         # print('{}: ensembl gpos {} != my gpos {}\n'.format(ind, ens_start, row['Sites']))
#     # else:
#     #     print('Okay\n')
# print('{} bad cDNA->gDNA conversions out of {}'.format(len(bad_indices), len(df_cheui_glori)))
#
# df_cheui_glori_filtered = df_cheui_glori.drop(bad_indices)
# df_cheui_glori_filtered.to_csv(cheui_file+'.glori.filtered', sep='\t', index=False)