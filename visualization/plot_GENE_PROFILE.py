import os
from Bio import SeqIO
import pysam
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
# mpl.rcParams['figure.dpi'] = 600
# mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['font.size'] = 5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
# fig_kwargs = dict(format=FMT, dpi=1200)
#######################################################################

def load_genome_reference(ref_file, chrs):
    print(f'Parsing genome reference {os.path.basename(ref_file)}...')
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        if record.id in chrs:
            ref[record.id] = str(record.seq)
    return ref

source_data_dir = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE'

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
bam_file = os.path.join(source_data_dir, 'chr6_31815543_31817946.mAFiA.reads.6motifs.bam')
img_out = '/home/adrian/NCOMMS_revision/images/GENE_PROFILE'

os.makedirs(img_out, exist_ok=True)
########################################################################################################################
### histogram of p(m6A) at single sites ################################################################################
########################################################################################################################
sel_chrom = '6'
sel_chromStart = 31815543
sel_chromEnd = 31817946

ref = load_genome_reference(ref_file, [sel_chrom])

ref_pos_mod_probs = {}
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for col in bam.pileup(sel_chrom, sel_chromStart, sel_chromEnd, truncate=True):
        if ref[sel_chrom][col.pos]=='A':
            # print("\ncoverage at base %s = %s" % (col.pos, col.n))
            all_mod_probs = []
            for pileupread in col.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # print('\tbase in read %s = %s' %
                    #       (pileupread.alignment.query_name,
                    #        pileupread.alignment.query_sequence[pileupread.query_position]))

                    for mod_pos, mod_probs in pileupread.alignment.modified_bases_forward.get(('N', 0, 21891), []):
                        if mod_pos==pileupread.query_position:
                            all_mod_probs.append(mod_probs/255.0)

            if len(all_mod_probs):
                ref_pos_mod_probs[col.pos] = all_mod_probs

## plot histogram ###
thresh_mod = 0.5
sel_pos = [31816312, 31816450, 31817452]
for ref_pos, mod_probs in ref_pos_mod_probs.items():
    # if len(mod_probs)>=400:
    if ref_pos in sel_pos:
        mod_ratio = int(np.mean(np.array(mod_probs)>=thresh_mod) * 100)
        plt.figure(figsize=(3*cm, 2*cm))
        plt.hist(mod_probs, bins=100, range=[0, 1])
        plt.xlim([-0.01, 1.01])
        plt.axvline(x=0.5, c='r', linestyle='--', alpha=0.5, linewidth=0.5)
        plt.title(f'chr{sel_chrom}: {ref_pos}\nS={mod_ratio}%')
        plt.savefig(os.path.join(img_out, f'hist_mod_probs_chr{sel_chrom}_{ref_pos}.{FMT}'), **fig_kwargs)
plt.close('all')

########################################################################################################################
### filter bam #########################################################################################################
########################################################################################################################
# def generate_mm_ml_tags(read_mods, mod_base='N', mod_code='21891'):
#     dists = [read_mods[0][0]] + list(np.diff([mod[0] for mod in read_mods]) - 1)
#     mod_probs = [mod[1] for mod in read_mods]
#     strand = '+'
#     mm_tag = mod_base + strand + mod_code + ',' + ','.join([str(d) for d in dists]) + ';'
#     ml_tag = mod_probs
#     return mm_tag, ml_tag
#
# sel_motifs = [
#     'AGACT',
#     'GAACT',
#     'GGACA',
#     'GGACC',
#     'GGACT',
#     'TGACT',
# ]
#
# out_bam_file = bam_file.replace('merged', f'chr{sel_chrom}_{sel_chromStart}_{sel_chromEnd}').replace('reads.bam', 'reads.6motifs.bam')
# with pysam.AlignmentFile(bam_file, 'rb') as bam_in:
#     with pysam.AlignmentFile(out_bam_file, 'wb', template=bam_in) as bam_out:
#         for read in bam_in.fetch(sel_chrom, sel_chromStart, sel_chromEnd):
#             old_read_mods = read.modified_bases_forward.get(('N', 0, 21891), [])
#             if len(old_read_mods):
#                 dict_aligned_pos = {read_pos: ref_pos for read_pos, ref_pos in read.get_aligned_pairs() if
#                                     ref_pos is not None}
#                 new_read_mods = []
#                 for mod_pos, mod_probs in old_read_mods:
#                     motif = ref[sel_chrom][dict_aligned_pos[mod_pos]-2:dict_aligned_pos[mod_pos]+3]
#                     if motif in sel_motifs:
#                         # print(motif)
#                         new_read_mods.append((mod_pos, mod_probs))
#                 if len(new_read_mods):
#                     new_mm, new_ml = generate_mm_ml_tags(new_read_mods)
#                     read.set_tag('MM', new_mm)
#                     read.set_tag('ML', new_ml)
#                     bam_out.write(read)