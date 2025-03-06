import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 10
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
# FMT = 'png'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
#######################################################################
import matplotlib.pyplot as plt
import os

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'

ds_x = 'HEK293_ctrl_RTA'
# ds_y = 'HEK293_TRUB1_kd_RTA'
# ds_y = 'HEK293_PUS1_kd_RTA'
ds_y = 'HEK293_PUS7_kd_RTA'

display_x = ds_x.lstrip('HEK293_').rstrip('_RTA')
display_y = ds_y.lstrip('HEK293_').rstrip('_RTA')

img_out = '/home/adrian/img_out/HEK293_psU_kds'
os.makedirs(img_out, exist_ok=True)

bedmethyl_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'frequency'
]

comp_mod_name = '17802'
dict_mod_display = {
    'a': 'm^6A',
    '17802': '\psi'
}

df_x = pd.read_csv(os.path.join(base_dir, ds_x, 'modkit042.cov10.bedmethyl'), sep='\t',
                   usecols=[0, 1, 2, 3, 4, 5, 10], names=bedmethyl_fields,
                   dtype={'chrom': str})
df_y = pd.read_csv(os.path.join(base_dir, ds_y, 'modkit042.cov10.bedmethyl'), sep='\t',
                   usecols=[0, 1, 2, 3, 4, 5, 10], names=bedmethyl_fields,
                   dtype={'chrom': str})

df_x_mod = df_x[df_x['name'] == comp_mod_name]
df_y_mod = df_y[df_y['name'] == comp_mod_name]
df_merged = pd.merge(df_x_mod, df_y_mod, how='inner', on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand'])
df_merged_filtered = df_merged[~((df_merged['frequency_x'] == 0.0) * (df_merged['frequency_y'] == 0.0))]
df_merged_filtered['delta'] = df_merged_filtered['frequency_y'] - df_merged_filtered['frequency_x']

num_sites = len(df_merged_filtered)

plt.figure(figsize=(10*cm, 5*cm))

mat_z, edges_x, edges_y = np.histogram2d(df_merged_filtered['frequency_y'], df_merged_filtered['frequency_x'],
                                         bins=20, range=[[0, 100], [0, 100]])
centers_x = 0.5 * (edges_x[1:] + edges_x[:-1])
centers_y = centers_x

plt.subplot(1, 2, 1)
xylim = [0, 100]
# plt.scatter(df_merged_filtered['frequency_x'], df_merged_filtered['frequency_y'], s=1)
plt.imshow(np.log10(mat_z+1), extent=xylim+xylim, origin='lower', vmin=0, vmax=3)
plt.plot([0, 100], [0, 100], 'r--')
plt.xticks(np.linspace(0, 100, 5))
plt.yticks(np.linspace(0, 100, 5))
plt.xlim(xylim)
plt.ylim(xylim)
plt.xlabel(display_x)
plt.ylabel(display_y)
# plt.savefig(os.path.join(img_out, f'scatter_{comp_mod_name}_{ds_x}_{ds_y}.png'), bbox_inches='tight')

plt.subplot(1, 2, 2)
plt.hist(df_merged_filtered['delta'], range=[-100, 100], bins=100)
plt.yscale('log')
plt.axvline(x=0, c='r', ls='--')
plt.xlabel('$\Delta$S(y-x)')
plt.ylabel('Site count')
# plt.title(f'{ds_y} - {ds_x}')

plt.suptitle(f'{num_sites} ${dict_mod_display[comp_mod_name]}$ sites')
plt.tight_layout()
plt.savefig(os.path.join(img_out, f'{display_x}_{display_y}_{comp_mod_name}.{FMT}'), **fig_kwargs)

### motif frequency ###
from Bio import SeqIO
from Bio.Seq import Seq
ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, "fasta"):
    ref[record.id] = record.seq

thresh_delta = -50
df_signif = df_merged_filtered[df_merged_filtered['delta'] < thresh_delta]
collected_motifs = []
for _, this_row in df_signif.iterrows():
    this_chrom, this_chromStart, this_chromEnd, this_strand = this_row[['chrom', 'chromStart', 'chromEnd', 'strand']]
    this_motif = ref[this_chrom][this_chromStart-2:this_chromStart+3]
    if this_strand == '-':
        this_motif = this_motif.reverse_complement()
    collected_motifs.append(str(this_motif))

import logomaker as lm
counts_mat = lm.alignment_to_matrix(collected_motifs)

plt.figure(figsize=(5*cm, 5*cm))
lm.Logo(counts_mat)
plt.xticks(np.arange(5), np.arange(-2, 3))
plt.xlabel('Position')
plt.ylabel('Count')
plt.title(f'{display_y}, $\Delta$S < {thresh_delta}')
plt.savefig(os.path.join(img_out, f'motifs_DeltaS{thresh_delta}_{display_x}_{display_y}_{comp_mod_name}.{FMT}'), **fig_kwargs)

plt.close('all')