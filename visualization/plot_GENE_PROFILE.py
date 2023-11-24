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
#######################################################################

def load_genome_reference(ref_file, chrs):
    print(f'Parsing genome reference {os.path.basename(ref_file)}...')
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        if record.id in chrs:
            ref[record.id] = str(record.seq)
    return ref

ref_file = '/home/adrian/Data/GRCh38_102/GRCh38_102.fa'
bam_file = '/home/adrian/NCOMMS_revision/source_data/GENE_PROFILE/merged.mAFiA.reads.bam'
img_out = '/home/adrian/NCOMMS_revision/images/GENE_PROFILE'

chrom = '6'
chromStart = 31815543
chromEnd = 31817946

ref = load_genome_reference(ref_file, [chrom])

ref_pos_mod_probs = {}
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for col in bam.pileup(chrom, chromStart, chromEnd, truncate=True):
        if ref[chrom][col.pos]=='A':
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

### plot histogram ###
thresh_mod = 0.5
# sel_pos = [31816450, 31817110, 31817566]
for ref_pos, mod_probs in ref_pos_mod_probs.items():
    if len(mod_probs)>=400:
    # if ref_pos in sel_pos:
        mod_ratio = int(np.mean(np.array(mod_probs)>=thresh_mod) * 100)
        plt.figure(figsize=(3*cm, 2*cm))
        plt.hist(mod_probs, bins=100, range=[0, 1])
        plt.xlim([-0.01, 1.01])
        plt.axvline(x=0.5, c='r', linestyle='--', alpha=0.5, linewidth=0.5)
        plt.title(f'chr{chrom}: {ref_pos}\nS={mod_ratio}%')
        plt.savefig(os.path.join(img_out, f'hist_mod_probs_chr{chrom}_{ref_pos}.{FMT}'), **fig_kwargs)
plt.close('all')