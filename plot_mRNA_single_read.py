import os
HOME = os.path.expanduser('~')
from glob import glob
import argparse
import pandas as pd
import numpy as np
import random
random.seed(0)
from random import sample
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

results_dir = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results'
train_dataset = 'WUE_batches1+2'
test_dataset = '100_WT_0_IVT'

img_out = os.path.join(HOME, 'img_out/MAFIA', os.path.basename('HEK293_single_read_{}'.format(train_dataset)))
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
    if counts>10:
        break
    this_site = unique_sites[i]
    next_site = unique_sites[i+1]
    if ((next_site-this_site)<=5):
        common_chr = set.intersection(set(df[df['Sites']==this_site]['Chr'].unique()), set(df[df['Sites']==next_site]['Chr'].unique()))
        if len(common_chr)>0:
            for this_common_chr in list(common_chr):
                print(this_common_chr, this_site, next_site)
                df_this_site = df[(df['Chr']==this_common_chr) * (df['Sites']==this_site)]
                df_next_site = df[(df['Chr']==this_common_chr) * (df['Sites']==next_site)]
                gene_name = df_this_site['Gene'].values[0]

                ### single-read view ###
                common_reads = list(set.intersection(set(df_this_site['read_id']), set(df_next_site['read_id'])))
                sample_reads = sample(common_reads, num_samples)

                mat_mod_prob = np.zeros([num_samples, extent])
                ids = []
                for i, read in enumerate(sample_reads):
                    ids.append(read.split('-')[-1])
                    mat_mod_prob[i, min_pos] = df_this_site[df_this_site['read_id']==read]['mod_prob'].values[0]
                    mat_mod_prob[i, min_pos+next_site-this_site] = df_next_site[df_next_site['read_id']==read]['mod_prob'].values[0]

                xtick_pos = [min_pos, min_pos+next_site-this_site]
                xticks = [this_site, next_site]

                fig_single_read = plt.figure(figsize=(5, 10))
                ax = fig_single_read.add_subplot(1, 1, 1)
                im = ax.imshow(mat_mod_prob, vmin=vmin, vmax=1)
                ax.set_xticks(xtick_pos, xticks, fontsize=8)
                ax.set_yticks(np.arange(num_samples), ids, fontsize=8)
                ax.set_xlabel('Aligned pos (NTs)', fontsize=12)
                ax.set_ylabel('Read IDs', fontsize=12)
                ax.set_title('{}: {}-{}\n{}'.format(this_common_chr, this_site-min_pos, this_site-min_pos+extent, gene_name), fontsize=15)

                fig_single_read.tight_layout(rect=[0.1, 0.2, 0.8, 0.8])
                fig_single_read.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)
                cax = fig_single_read.add_axes([0.75, 0.3, 0.05, 0.4])
                plt.colorbar(im, cax=cax)
                plt.ylabel('P(m6A)', fontsize=12, rotation=-90, labelpad=30)
                plt.yticks(cax_yticks, cax_yticks)

                fig_single_read.savefig(os.path.join(img_out, 'single_read_{}_{}_{}_{}.png'.format(this_common_chr, this_site, next_site, gene_name)), bbox_inches='tight')

                counts += 1

                ### prob. distribution ###
                this_site_mod_probs = df_this_site['mod_prob'].values
                next_site_mod_probs = df_next_site['mod_prob'].values

                this_site_GLORI = df_this_site['Ratio'].values[0]
                next_site_GLORI = df_next_site['Ratio'].values[0]

                this_site_mod_ratio = np.mean(this_site_mod_probs>=thresh_mod_prob)
                next_site_mod_ratio = np.mean(next_site_mod_probs>=thresh_mod_prob)

                plt.figure(figsize=(10, 6))
                plt.subplot(1, 2, 1)
                plt.hist(this_site_mod_probs, bins=50, range=[0, 1])
                plt.ylabel('Counts', fontsize=12)
                plt.xlabel('Single-NT mod. prob.', fontsize=12)
                plt.title('{}: {}\nPred. mod. ratio {:.2f}\nGLORI {:.2f}'.format(this_common_chr, this_site, this_site_mod_ratio, this_site_GLORI), fontsize=15)
                plt.subplot(1, 2, 2)
                plt.hist(next_site_mod_probs, bins=50, range=[0, 1])
                plt.ylabel('Counts', fontsize=12)
                plt.xlabel('Single-NT mod. prob.', fontsize=12)
                plt.title('{}: {}\nPred. mod. ratio {:.2f}\nGLORI {:.2f}'.format(this_common_chr, next_site, next_site_mod_ratio, next_site_GLORI), fontsize=15)

                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                plt.suptitle(gene_name, fontsize=15)

                plt.savefig(os.path.join(img_out, 'hist_{}_{}_{}_{}.png'.format(this_common_chr, this_site, next_site, gene_name)), bbox_inches='tight')