import pandas as pd
import pysam
from tqdm import tqdm
import re
import os
import numpy as np

ds = 'Diet'
cond = 'M3KO_WD'

dict_mod_code = {
    'm6A': 21891,
    'psi': 17802,
}

top_mods = 5

def get_single_read_mod_probs(in_read, num_top_mods=5):
    dict_query_to_ref = {query_pos: ref_pos for query_pos, ref_pos in
                         this_read.get_aligned_pairs(matches_only=True)}
    mod_top_refPos_modProb = {}
    for this_mod in dict_mod_code.keys():
        mod_bases = in_read.modified_bases.get(('N', 0, dict_mod_code[this_mod]), [])
        if len(mod_bases):
            top_mod_bases = [mod_bases[this_ind]
                             for this_ind in np.argsort([tup[1] for tup in mod_bases])[-num_top_mods:][::-1]]
            mod_top_refPos_modProb[this_mod] = [(dict_query_to_ref[this_tup[0]], round(this_tup[1]/255.0, 3))
                                                for this_tup in top_mod_bases]
    return mod_top_refPos_modProb


polyA_dir = f"/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/polyA"
polyA_file = os.path.join(polyA_dir, f"{cond}_read2polyA_length.txt")
if ds == 'TAC':
    bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/{cond}/chrALL.mAFiA.reads.bam'
else:
    bam_file = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/{ds}/{cond}_merged/chrALL.mAFiA.reads.bam'
gene_bed = '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/gene.GRCm38.102.bed'

out_dir_polyA = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA'
os.makedirs(out_dir_polyA, exist_ok=True)

out_dir_top_mod = f'/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/top_mods'
os.makedirs(out_dir_top_mod, exist_ok=True)

bed_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'source',
    'biotype',
    'frame',
    'description'
]

df_polyA = pd.read_csv(polyA_file, sep='\t', names=['read_id', 'polyA_length'])
read_ids_with_polyA = list(df_polyA['read_id'].values)

df_gene = pd.read_csv(gene_bed, sep='\t', names=bed_fields)
df_gene = df_gene[df_gene['source'] == 'ensembl_havana']
df_gene['gene'] = [re.search('gene_name "(.*?)";', desc).group(1) for desc in df_gene['description']]
gene_to_chrom = {this_row['gene']: this_row['chrom'] for _, this_row in df_gene.iterrows()}

gene_to_read = {}
read_mod_probs = {}
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for _, this_row in tqdm(df_gene.iterrows()):
        if this_row['strand'] == '+':
            flag_required = 0
        else:
            flag_required = 16
        collected_read_ids = []
        for this_read in bam.fetch(this_row['chrom'], this_row['chromStart'], this_row['chromEnd']):
            if this_read.flag == flag_required:
                collected_read_ids.append(this_read.query_name)
                read_mod_probs[this_read.query_name] = get_single_read_mod_probs(this_read, num_top_mods=top_mods)
        if len(collected_read_ids):
            gene_to_read[this_row['gene']] = collected_read_ids

read_to_gene = {
    read: gene
    for gene, reads in gene_to_read.items()
    for read in reads
}

df_polyA['gene'] = [read_to_gene.get(this_read) for this_read in df_polyA['read_id']]
df_polyA_sel = df_polyA[~df_polyA['gene'].isna()]

df_polyA_sel.to_csv(os.path.join(out_dir_polyA, f'polyA_reads_annotated_{ds}_{cond}.tsv'), sep='\t', index=False)
print(f'{len(df_polyA_sel)} reads written with polyA length')

top_mod_fields = [
    'read_id',
    'gene',
    'chrom',
    'refPos_modProb_m6A',
    'refPos_modProb_psi'
]

df_top_mod_probs = pd.DataFrame(
    {
        'read_id': read_mod_probs.keys(),
        'refPos_modProb_m6A': [this_val.get('m6A') for this_val in read_mod_probs.values()],
        'refPos_modProb_psi': [this_val.get('psi') for this_val in read_mod_probs.values()],
    }
)
df_top_mod_probs['gene'] = [read_to_gene.get(this_read_id) for this_read_id in df_top_mod_probs['read_id']]
df_top_mod_probs['chrom'] = [gene_to_chrom.get(this_gene) for this_gene in df_top_mod_probs['gene']]
df_top_mod_probs[top_mod_fields].to_csv(os.path.join(out_dir_top_mod, f'top_{top_mods}_mods_per_read_{ds}_{cond}.tsv'),
                                        sep='\t', index=False, float_format='%.3f')

print(f'{len(df_top_mod_probs)} reads written with top mod. locations')
