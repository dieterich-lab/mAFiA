import os
import pysam
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.stats import kstest

# chrom, chromStart, chromEnd, name, strand = ('1', 4776649, 4776650,	'm6A', '-')
chrom, chromStart, chromEnd, name, strand = ('11', 95680478, 95680479, 'm6A', '+')

ds = 'HFpEF'
conditions = ['ctrl_merged', 'HFpEF_merged']

cond_names = {
    this_cond: this_cond.rstrip('_merged') for this_cond in conditions
}
cond_colors = {
    this_cond: this_color for this_cond, this_color in zip(conditions, ['b', 'r'])
}

mod_code = 21891 if name=='m6A' else 17802

res_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart'

bam_files = {
    this_cond: os.path.join(res_dir, ds, this_cond, 'chrALL.mAFiA.reads.bam')
    for this_cond in conditions
}

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

mod_probs = {}
for this_cond in conditions:
    with pysam.AlignmentFile(bam_files[this_cond], 'rb') as bam:
        mod_probs[this_cond] = collect_mod_prob(bam, chrom, chromStart, chromEnd, strand)

ks_stat, ks_pval = kstest(*list(mod_probs.values()))


xylim = [-0.01, 1.01]

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
for this_cond in conditions:
    this_cond_mod_probs = mod_probs[this_cond]
    modRatio, coverage, confidence = calc_modRatio_coverage_confidence(this_cond_mod_probs)
    plt.hist(this_cond_mod_probs, bins=50, range=[0, 1],
             fc=cond_colors[this_cond], alpha=0.5,
             label=f'{cond_names[this_cond]}\nN={coverage}, S={modRatio}%, conf={confidence}')
plt.xlim(xylim)
plt.xlabel(f'$P({{{dict_mod_display[name]}}})$', fontsize=12)
plt.ylabel('Num. reads', fontsize=12)
# plt.title(rf"$S_{'SHAM'}$={modRatio_sham}, $S_{'TAC'}$={modRatio_tac}")
plt.legend()
plt.subplot(1, 2, 2)
for this_cond in conditions:
    plt.ecdf(mod_probs[this_cond], c=cond_colors[this_cond], label=cond_names[this_cond])
plt.title(f'KS stat. {ks_stat:.3f}, p-val {ks_pval:.3f}')
plt.xlim(xylim)
plt.ylim(xylim)
plt.xlabel(f'$P({{{dict_mod_display[name]}}})$', fontsize=12)
plt.ylabel('cdf', fontsize=12)
plt.legend(loc='lower right')
plt.suptitle(f'chr{chrom}: {chromEnd}, {strand}', fontsize=15)
plt.savefig(os.path.join(img_out, f'chr{chrom}_{chromEnd}_strand{strand}.png'), bbox_inches='tight')