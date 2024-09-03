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
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['font.family'] = 'Arial'
FMT = 'pdf'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=1200)
#######################################################################

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def load_genome_reference(ref_file, chrs=None):
    if chrs is None:
        chrs = [str(chr) for chr in range(1, 23)] + ['X']
    print(f'Parsing genome reference {os.path.basename(ref_file)}...')
    ref = {}
    for record in SeqIO.parse(ref_file, 'fasta'):
        if record.id in chrs:
            ref[record.id] = str(record.seq)
    return ref

ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
bam_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/HEK293/WT_P2/chrALL.mAFiA.reads.bam'

img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure1'
os.makedirs(img_out, exist_ok=True)

ref = load_genome_reference(ref_file)

# sel_chrom = '9'
# sel_chromStart = 137108465
# sel_chromEnd = 137108663

sel_chrom = '1'
sel_chromStart = 151008481
sel_chromEnd = 151008658

ref_pos_mod_probs = {}
with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for col in bam.pileup(sel_chrom, sel_chromStart, sel_chromEnd, truncate=True):
        if ref[sel_chrom][col.pos]=='A':
            mod_code = 21891
        elif ref[sel_chrom][col.pos]=='T':
            mod_code = 17802
        else:
            continue

        all_mod_probs = []
        for pileupread in col.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                for mod_pos, mod_probs in pileupread.alignment.modified_bases_forward.get(('N', 0, mod_code), []):
                    if mod_pos==pileupread.query_position:
                        all_mod_probs.append(mod_probs/255.0)

            if len(all_mod_probs):
                ref_pos_mod_probs[col.pos] = all_mod_probs

## plot histogram ###
thresh_mod = 0.5
sel_pos = [151008502, 151008514, 151008569, 151008639]
for ref_pos, mod_probs in ref_pos_mod_probs.items():
    if ref_pos in sel_pos:
        mod_ratio = int(np.mean(np.array(mod_probs)>=thresh_mod) * 100)
        ref5mer = ref[sel_chrom][ref_pos-2:ref_pos+3]
        if ref5mer[2]=='A':
            mod = 'm6A'
            color = 'red'
        elif ref5mer[2]=='T':
            mod = 'psi'
            color = 'blue'
        plt.figure(figsize=(2*cm, 2*cm))
        counts, bins, _ = plt.hist(mod_probs, bins=20, range=[0, 1], facecolor=color)
        plt.xlim([-0.01, 1.01])
        plt.xticks([0, 0.5, 1])
        # plt.xlabel(f'$P({dict_mod_display[mod]})$')
        # plt.ylabel('Count NTs')
        plt.axvline(x=0.5, c='r', linestyle='--', alpha=0.5, linewidth=0.5)
        # plt.title(f"chr{sel_chrom}: {ref_pos+1}\n{ref5mer.replace('T', 'U')}, S={mod_ratio}%")
        plt.title(rf"$S_{{{dict_mod_display[mod]}}}$ = {mod_ratio}%")
        plt.savefig(os.path.join(img_out, f'hist_mod_probs_chr{sel_chrom}_{ref_pos}_{ref5mer}.{FMT}'), **fig_kwargs)
plt.close('all')