import os
import pysam
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.stats import kstest

# chrom = '1'
# chromStart = 37417732
# chromEnd = 37417733
# name = 'm6A'
# strand = '-'

# chrom = '1'
# chromStart = 37417500
# chromEnd = 37417501
# name = 'psi'
# strand = '-'

# chrom = '6'
# chromStart = 121223109
# chromEnd = 121223110
# name = 'm6A'
# strand = '+'

# chrom = '6'
# chromStart = 121223106
# chromEnd = 121223107
# name = 'psi'
# strand = '+'

chrom, chromStart, chromEnd, name, strand = ('1', 4776649, 4776650,	'm6A', '-')

mod_code = 21891 if name=='m6A' else 17802

bam_sham_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_day21/chrALL.mAFiA.reads.bam'
bam_tac_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_day21/chrALL.mAFiA.reads.bam'

img_out = '/home/adrian/img_out/single_site_ks_test'
os.makedirs(img_out, exist_ok=True)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

def collect_mod_prob(in_bam, in_chrom, in_chromStart, in_chromEnd, in_strand):
    out_mod_probs = []
    flag_require = 0 if in_strand=='+' else 16
    for pileupcolumn in in_bam.pileup(in_chrom, in_chromStart, in_chromEnd, truncate=True, flag_require=flag_require):
        if pileupcolumn.pos==in_chromStart:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    this_read_mod_probs = [score for (pos, score) in pileupread.alignment.modified_bases.get(('N', 0, mod_code), [])
                     if pos == pileupread.query_position]
                    out_mod_probs.extend(this_read_mod_probs)
    return np.array(out_mod_probs) / 255.0


def calc_modRatio_coverage_confidence(in_mod_probs):
    out_modRatio = np.mean(in_mod_probs >= 0.5) * 100.0
    out_coverage = len(in_mod_probs)
    out_confidence = ((in_mod_probs < 0.25).sum() + (in_mod_probs >= 0.75).sum()) / out_coverage * 100.0
    return np.round(out_modRatio, 1), np.round(out_coverage, 1), np.round(out_confidence, 1)


with pysam.AlignmentFile(bam_sham_file, 'rb') as bam_sham:
    with pysam.AlignmentFile(bam_tac_file, 'rb') as bam_tac:
        mod_probs_sham = collect_mod_prob(bam_sham, chrom, chromStart, chromEnd, strand)
        mod_probs_tac = collect_mod_prob(bam_tac, chrom, chromStart, chromEnd, strand)

# mod_probs_sham = np.array(mod_probs_sham) / 255.0
# mod_probs_tac = np.array(mod_probs_tac) / 255.0
ks_stat, ks_pval = kstest(mod_probs_sham, mod_probs_tac)
modRatio_sham, coverage_sham, confidence_sham = calc_modRatio_coverage_confidence(mod_probs_sham)
modRatio_tac, coverage_tac, confidence_tac = calc_modRatio_coverage_confidence(mod_probs_tac)

xlim = [-0.01, 1.01]

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(mod_probs_sham, bins=50, range=[0, 1], fc='b', alpha=0.5, label=f'SHAM, N={coverage_sham}, S={modRatio_sham}%, conf={confidence_sham}')
plt.hist(mod_probs_tac, bins=50, range=[0, 1], fc='r', alpha=0.5, label=f'TAC, N={coverage_tac}, S={modRatio_tac}%, conf={confidence_tac}')
plt.xlim(xlim)
plt.xlabel(f'$P({{{dict_mod_display[name]}}})$', fontsize=12)
plt.ylabel('Num. reads', fontsize=12)
# plt.title(rf"$S_{'SHAM'}$={modRatio_sham}, $S_{'TAC'}$={modRatio_tac}")
plt.legend()
plt.subplot(1, 2, 2)
plt.ecdf(mod_probs_sham, c='b', label='SHAM')
plt.ecdf(mod_probs_tac, c='r', label='TAC')
plt.title(f'KS stat. {ks_stat:.3f}, p-val {ks_pval:.3f}')
plt.xlim(xlim)
plt.xlabel(f'$P({{{dict_mod_display[name]}}})$', fontsize=12)
plt.ylabel('cdf', fontsize=12)
plt.legend(loc='upper left')
plt.suptitle(f'chr{chrom}: {chromEnd}, {strand}', fontsize=15)
# plt.savefig(os.path.join(img_out, f'chr{chrom}_{chromEnd}_strand{strand}.png'), bbox_inches='tight')