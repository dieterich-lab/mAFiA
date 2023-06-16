import os
HOME = os.path.expanduser('~')
import pandas as pd
import numpy as np
import random
random.seed(0)
from random import sample

#######################################################################
import matplotlib as mpl
import matplotlib.pyplot as plt

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

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
train_dataset = 'ISA-WUE'
test_dataset = '100_WT_0_IVT'

img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('mRNA_train_{}_test_HEK293_single_read'.format(train_dataset)))
if not os.path.exists(img_out):
    os.makedirs(img_out, exist_ok=True)

path = os.path.join(results_dir, 'res_train_{}_test_{}.tsv.merged'.format(train_dataset, test_dataset))
df = pd.read_csv(path, sep='\t').rename(columns={'Unnamed: 0': 'index'})

num_samples = 50
min_pos = 2
counts = 0
vmin = 0.8
cax_yticks = np.arange(vmin, 1.01, 0.1)
extent = 10

thresh_mod_prob = 0.5

### choose nearby sites ###
unique_sites = df['Sites'].unique()
unique_sites.sort()
for i in range(len(unique_sites)-1):
    # if counts>10:
    #     break
    this_site = unique_sites[i]
    next_site = unique_sites[i+1]
    if ((next_site-this_site)<=5):
        common_chr = set.intersection(set(df[df['Sites']==this_site]['Chr'].unique()), set(df[df['Sites']==next_site]['Chr'].unique()))

        if (len(common_chr)>0):
            for this_common_chr in list(common_chr):
                print(this_common_chr, this_site, next_site)
                df_this_site = df[(df['Chr']==this_common_chr) * (df['Sites']==this_site)]
                df_next_site = df[(df['Chr']==this_common_chr) * (df['Sites']==next_site)]
                common_gene = set.intersection(set(df_this_site['Gene'].unique()),
                                               set(df_next_site['Gene'].unique()))
                if len(common_gene)>1:
                    continue
                gene_name = df_this_site['Gene'].values[0]

                ### single-read view ###
                common_reads = list(set.intersection(set(df_this_site['read_id']), set(df_next_site['read_id'])))
                if len(common_reads)<num_samples:
                    continue
                sample_reads = sample(common_reads, num_samples)

                mat_mod_prob = np.zeros([num_samples, extent])
                ids = []
                for i, read in enumerate(sample_reads):
                    ids.append(read.split('-')[-1])
                    mat_mod_prob[i, min_pos] = df_this_site[df_this_site['read_id']==read]['mod_prob'].values[0]
                    mat_mod_prob[i, min_pos+next_site-this_site] = df_next_site[df_next_site['read_id']==read]['mod_prob'].values[0]

                xtick_pos = [min_pos, min_pos+next_site-this_site]
                xticks = [this_site, next_site]

                ytick_pos = np.concatenate([[0], np.arange(9, num_samples, 10)])
                yticks = [f'read {i+1}' for i in ytick_pos]


                fig_single_read = plt.figure(figsize=(8*cm, 10*cm))
                ax = fig_single_read.add_subplot(1, 1, 1)
                im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1)
                ax.set_xticks(xtick_pos, xticks, rotation=-90)
                ax.set_yticks(ytick_pos, yticks)
                # ax.set_xlabel('Genome Pos.')
                # ax.set_ylabel('Read IDs')
                ax.set_title(f'{gene_name}\n{this_common_chr}: {this_site-min_pos}-{this_site-min_pos+extent}')

                # fig_single_read.tight_layout(rect=[0.1, 0.2, 0.8, 0.8])
                # fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
                cax = fig_single_read.add_axes([0.7, 0.3, 0.03, 0.4])
                plt.colorbar(im, cax=cax)
                # plt.ylabel('P(m6A)', rotation=-90)
                plt.yticks(cax_yticks, cax_yticks)

                fig_single_read.savefig(os.path.join(img_out, f'single_read_{this_common_chr}_{this_site}_{next_site}_{gene_name}.{FMT}'), **fig_kwargs)

                counts += 1

                ### prob. distribution ###
                this_site_mod_probs = df_this_site['mod_prob'].values
                next_site_mod_probs = df_next_site['mod_prob'].values

                this_site_GLORI = df_this_site['Ratio'].values[0]
                next_site_GLORI = df_next_site['Ratio'].values[0]

                this_site_mod_ratio = np.mean(this_site_mod_probs>=thresh_mod_prob)
                next_site_mod_ratio = np.mean(next_site_mod_probs>=thresh_mod_prob)

                this_site_counts, bin_edges = np.histogram(this_site_mod_probs, bins=50, range=[0, 1])
                next_site_counts, _ = np.histogram(next_site_mod_probs, bins=50, range=[0, 1])

                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                xticks = np.arange(0, 1.01, 0.25)

                plt.figure(figsize=(8*cm, 10*cm))
                plt.subplot(2, 1, 1)
                plt.bar(bin_centers, this_site_counts, width=0.02, label=f'mAFiA {this_site_mod_ratio:.2f}\nGLORI {this_site_GLORI:.2f}')
                ymax = (np.max(this_site_counts) // 10 + 2) * 10
                yticks = np.arange(0, ymax+1, 10)
                plt.xlim([-0.01, 1.01])
                plt.ylabel('Read Counts')
                # plt.xticks(xticks)
                plt.xticks([])
                plt.yticks(yticks)
                # plt.xlabel('Single-NT mod. prob.')
                plt.axvline(x=0.5, c='r', alpha=0.5, linewidth=0.5, linestyle='--')
                plt.legend(title=f'{this_common_chr}: {this_site}', loc='upper left', markerscale=0.01)
                # plt.title(f'{this_common_chr}: {this_site}')
                plt.subplot(2, 1, 2)
                plt.bar(bin_centers, next_site_counts, width=0.02, label=f'mAFiA {next_site_mod_ratio:.2f}\nGLORI {next_site_GLORI:.2f}')
                ymax = (np.max(next_site_counts) // 10 + 2) * 10
                yticks = np.arange(0, ymax+1, 10)
                plt.xlim([-0.01, 1.01])
                plt.ylabel('Read Counts')
                plt.xlabel('Modification probability')
                plt.xticks(xticks)
                plt.yticks(yticks)
                plt.axvline(x=0.5, c='r', alpha=0.5, linewidth=0.5, linestyle='--')
                plt.legend(title=f'{this_common_chr}: {next_site}', loc='upper left', markerscale=0.01)
                # plt.title(f'{this_common_chr}: {next_site}')

                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                plt.suptitle(gene_name)
                plt.savefig(os.path.join(img_out, f'hist_{this_common_chr}_{this_site}_{next_site}_{gene_name}.{FMT}'), **fig_kwargs)
                plt.close('all')