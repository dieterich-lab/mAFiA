from Bio import SeqIO
from tqdm import tqdm
import logomaker as lm
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
mpl.rcParams['font.size'] = 8
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
#######################################################################
import matplotlib.pyplot as plt
import os


def get_central_motif(in_df, span=2):
    all_motifs = []
    for _, this_row in tqdm(in_df.iterrows()):
        this_chrom, this_chromStart, this_chromEnd, this_strand = this_row[
            ['chrom', 'chromStart', 'chromEnd', 'strand']
        ]
        this_motif = ref[this_chrom][(this_chromStart-span):(this_chromStart+span+1)]
        if this_strand == '-':
            this_motif = this_motif.reverse_complement()
        all_motifs.append(str(this_motif))
    in_df['ref_motif'] = all_motifs
    return in_df


def plot_5mer_motif(in_mod_name, in_df, op, thresh_delta):
    ### motif frequency ###
    if op == '>=':
        df_signif = in_df[in_df['delta'].ge(thresh_delta)]
    elif op == '<':
        df_signif = in_df[in_df['delta'].lt(thresh_delta)]
    collected_motifs = []
    for _, this_row in df_signif.iterrows():
        this_chrom, this_chromStart, this_chromEnd, this_strand = this_row[
            ['chrom', 'chromStart', 'chromEnd', 'strand']]
        this_motif = ref[this_chrom][this_chromStart - 2:this_chromStart + 3]
        if this_strand == '-':
            this_motif = this_motif.reverse_complement()
        collected_motifs.append(str(this_motif))
    counts_mat = lm.alignment_to_matrix(collected_motifs)

    this_fig = plt.figure(figsize=(5 * cm, 5 * cm))
    lm.Logo(counts_mat)
    plt.xticks(np.arange(5), np.arange(-2, 3))
    plt.xlabel('Position')
    plt.ylabel('Count')
    plt.title(f'{display_y}, $\Delta$S(${dict_mod_display[in_mod_name]}$) {op} {thresh_delta}')
    plt.savefig(os.path.join(img_out, f'motifs_DeltaS{in_mod_name}{op}{thresh_delta}_{display_x}_{display_y}.{FMT}'),
                **fig_kwargs)
    this_fig.clf()


chemistry = 'RNA004'

### RNA004 ###
if chemistry == 'RNA004':
    base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'
    ds_x = 'HEK293_ctrl_RTA'

    sel_m6A_motifs = [
        'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
        'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
        'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
    ]

    ds_y = 'HEK293_TRUB1_kd_RTA'
    sel_psi_motifs = ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT']

    # ds_y = 'HEK293_PUS1_kd_RTA'
    # sel_psi_motifs = ['GTG', 'GTA', 'ATA', 'ATG']

    # ds_y = 'HEK293_PUS7_kd_RTA'
    # sel_psi_motifs = ['TGTAG']

    display_x = ds_x.lstrip('HEK293_').rstrip('_RTA')
    display_y = ds_y.lstrip('HEK293_').rstrip('_RTA')

### RNA002 ###
elif chemistry == 'RNA002':
    base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'
    ds_x = 'HEK293/WT_P2'
    ds_y = 'HEK293_TRUB1_OE/merged'

    display_x = 'WT'
    display_y = 'TRUB1_OE'

img_out = '/home/adrian/img_out/HEK293_psU_kds'
os.makedirs(img_out, exist_ok=True)

ref_file = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
ref = {}
for record in SeqIO.parse(ref_file, "fasta"):
    ref[record.id] = record.seq

if chemistry == 'RNA004':
    bed_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'score',
        'strand',
        'frequency'
    ]

    mod_names = ['17802', 'a']
    dict_mod_display = {
        'a': 'm^6A',
        '17802': '\psi'
    }

    df_x = pd.read_csv(os.path.join(base_dir, ds_x, 'modkit042.cov10.bedmethyl'), sep='\t',
                       usecols=[0, 1, 2, 3, 4, 5, 10], names=bed_fields,
                       dtype={'chrom': str})
    df_y = pd.read_csv(os.path.join(base_dir, ds_y, 'modkit042.cov10.bedmethyl'), sep='\t',
                       usecols=[0, 1, 2, 3, 4, 5, 10], names=bed_fields,
                       dtype={'chrom': str})

elif chemistry == 'RNA002':
    bed_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'score',
        'strand',
        'ref5mer',
        'coverage',
        'modRatio',
        'confidence'
    ]

    mod_names = ['psi', 'm6A']
    dict_mod_display = {
        'm6A': 'm^6A',
        'psi': '\psi'
    }

    thresh_cov = 10
    thresh_conf = 80.0

    df_x = pd.read_csv(os.path.join(base_dir, ds_x, 'chrALL.mAFiA.sites.bed'),
                       sep='\t', dtype={'chrom': str})
    df_x.rename(columns={'modRatio': 'frequency'}, inplace=True)
    df_x = df_x[
        (df_x['coverage'] >= thresh_cov)
        * (df_x['confidence'] >= thresh_conf)
    ]
    df_y = pd.read_csv(os.path.join(base_dir, ds_y, 'chrALL.mAFiA.sites.bed'),
                       sep='\t', dtype={'chrom': str})
    df_y = df_y[
        (df_y['coverage'] >= thresh_cov)
        * (df_y['confidence'] >= thresh_conf)
    ]
    df_y.rename(columns={'modRatio': 'frequency'}, inplace=True)


mod_df_merged_motif_filtered = {}
for mod_ind, mod_name in enumerate(mod_names):
    df_x_mod = df_x[df_x['name'] == mod_name]
    df_y_mod = df_y[df_y['name'] == mod_name]
    df_merged = pd.merge(df_x_mod, df_y_mod, how='inner', on=['chrom', 'chromStart', 'chromEnd', 'name', 'strand'])
    df_merged_filtered = df_merged[~((df_merged['frequency_x'] == 0.0) * (df_merged['frequency_y'] == 0.0))]
    df_merged_filtered['delta'] = df_merged_filtered['frequency_y'] - df_merged_filtered['frequency_x']

    if mod_name == '17802':
        if ds_y == 'HEK293_PUS1_kd_RTA':
            df_merged_filtered = get_central_motif(df_merged_filtered, span=1)
        else:
            df_merged_filtered = get_central_motif(df_merged_filtered, span=2)
        df_merged_motif_filtered = df_merged_filtered[df_merged_filtered['ref_motif'].isin(sel_psi_motifs)]
    elif mod_name == 'a':
        df_merged_filtered = get_central_motif(df_merged_filtered, span=2)
        df_merged_motif_filtered = df_merged_filtered[df_merged_filtered['ref_motif'].isin(sel_m6A_motifs)]

    mod_df_merged_motif_filtered[mod_name] = df_merged_motif_filtered


fig1 = plt.figure(figsize=(10 * cm, 10 * cm))
fig1_axes = fig1.subplots(2, 2)
for mod_ind, mod_name in enumerate(mod_names):
    df_merged_motif_filtered = mod_df_merged_motif_filtered[mod_name]
    num_sites = len(df_merged_motif_filtered)
    mat_z, edges_x, edges_y = np.histogram2d(df_merged_motif_filtered['frequency_y'], df_merged_motif_filtered['frequency_x'],
                                             bins=20, range=[[0, 100], [0, 100]])
    centers_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    centers_y = centers_x

    xylim = [0, 100]
    if mod_name == '17802':
        fig1_axes[mod_ind, 0].scatter(df_merged_motif_filtered['frequency_x'],
                                      df_merged_motif_filtered['frequency_y'],
                                      s=1)
    else:
        fig1_axes[mod_ind, 0].imshow(np.log10(mat_z+1), extent=xylim+xylim, origin='lower', vmin=0, vmax=3)
    fig1_axes[mod_ind, 0].plot([0, 100], [0, 100], 'r--')
    fig1_axes[mod_ind, 0].set_xticks(np.linspace(0, 100, 5))
    fig1_axes[mod_ind, 0].set_yticks(np.linspace(0, 100, 5))
    fig1_axes[mod_ind, 0].set_xlim(xylim)
    fig1_axes[mod_ind, 0].set_ylim(xylim)
    fig1_axes[mod_ind, 0].set_xlabel(display_x)
    fig1_axes[mod_ind, 0].set_ylabel(display_y)
    # plt.savefig(os.path.join(img_out, f'scatter_{comp_mod_name}_{ds_x}_{ds_y}.png'), bbox_inches='tight')
    fig1_axes[mod_ind, 0].set_title(f'{num_sites} ${dict_mod_display[mod_name]}$ sites')

    vec_delta = df_merged_motif_filtered['delta']
    vec_delta_pos = vec_delta[vec_delta >= 0]
    vec_delta_neg = vec_delta[vec_delta < 0]

    fig1_axes[mod_ind, 1].hist(vec_delta_pos, range=[0, 100], bins=20,
             density=True, histtype='step', edgecolor='r', alpha=0.75, label='+ve')
    fig1_axes[mod_ind, 1].hist(-vec_delta_neg, range=[0, 100], bins=20,
             density=True, histtype='step', edgecolor='b', alpha=0.75, label='-ve')
    fig1_axes[mod_ind, 1].set_yscale('log')
    fig1_axes[mod_ind, 1].set_xlabel(f'$\Delta$S(${dict_mod_display[mod_name]}$)')
    fig1_axes[mod_ind, 1].set_ylabel('Density')
    fig1_axes[mod_ind, 1].legend(loc='upper right')
    # plt.title(f'{ds_y} - {ds_x}')

    # if display_y.split('_')[1] == 'OE':
    #     if mod_name == 'psi':
    #         plot_5mer_motif(mod_name, df_merged_filtered, '>=', 25)
    #     elif mod_name == 'm6A':
    #         plot_5mer_motif(mod_name, df_merged_filtered, '<', -25)
    # elif display_y.split('_')[1] == 'kd':
    #     if mod_name == '17802':
    #         plot_5mer_motif(mod_name, df_merged_filtered, '<', -50)
    #     elif mod_name == 'a':
    #         plot_5mer_motif(mod_name, df_merged_filtered, '>=', 50)
fig1.tight_layout()
fig1.savefig(os.path.join(img_out, f'{display_x}_{display_y}.{FMT}'), **fig_kwargs)
plt.close('all')
