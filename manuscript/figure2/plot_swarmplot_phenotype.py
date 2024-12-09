import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 6
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
# FMT = 'png'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
#######################################################################
import matplotlib.pyplot as plt
import seaborn as sns


phenotype_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/mouse_phenotype_TAC_HFpEF.tsv'
df_phenotype = pd.read_csv(phenotype_file, sep='\t')
img_out = '/home/adrian/img_out/manuscript_psico_mAFiA/figure2'

HFrEF_conditions = ['Sham', 'TAC']
HFpEF_conditions = ['Ctrl', 'HFpEF']

df_HFrEF = df_phenotype[df_phenotype['Condition'].isin(HFrEF_conditions)]
df_HFpEF = df_phenotype[df_phenotype['Condition'].isin(HFpEF_conditions)]
df_condition = {
    'HFrEF': df_HFrEF,
    'HFpEF': df_HFpEF
}

cond_palette = {
    'HFrEF': ['gray', 'green'],
    'HFpEF': ['gray', 'red']
}

plt.figure(figsize=(6*cm, 6*cm))
for cond_ind, (this_cond, this_df_cond) in enumerate(df_condition.items()):

    plt.subplot(2, 2, 2*cond_ind+1)
    sns.stripplot(data=this_df_cond, x='Condition', y='Ejection Fraction',
                  hue='Condition', palette=cond_palette[this_cond], size=2)
    sns.barplot(data=this_df_cond, x='Condition', y='Ejection Fraction',
                  hue='Condition', palette=cond_palette[this_cond], width=0.5, alpha=0.5)
    plt.ylim([0, 100])
    plt.yticks([0, 50, 100])
    plt.xlabel(None)
    plt.ylabel('Ejection Fraction (%)')
    # if cond_ind == 1:
    #     plt.xlabel('Ejection Fraction (%)')
    # else:
    #     plt.xlabel(None)

    plt.subplot(2, 2, 2*cond_ind+2)
    sns.stripplot(data=this_df_cond, x='Condition', y='HW.BW (mg/g)',
                  hue='Condition', palette=cond_palette[this_cond], size=2)
    sns.barplot(data=this_df_cond, x='Condition', y='HW.BW (mg/g)',
                  hue='Condition', palette=cond_palette[this_cond], width=0.5, alpha=0.5)
    plt.ylim([3, 8])
    # plt.xticks(range(2), [])
    plt.xlabel(None)
    plt.ylabel('HW/BW (mg/g)')

    # if cond_ind == 1:
    #     plt.xlabel('HW/BW (mg/g)')
    # else:
    #     plt.xlabel(None)
plt.savefig(os.path.join(img_out, f'mouse_phenotype.{FMT}'), **fig_kwargs)

plt.figure(figsize=(4*cm, 4*cm))
this_cond = 'HFpEF'
sns.stripplot(data=this_df_cond, x='Condition', y='E/E`',
              hue='Condition', palette=cond_palette[this_cond], size=2)
sns.barplot(data=this_df_cond, x='Condition', y='E/E`',
              hue='Condition', palette=cond_palette[this_cond], width=0.5, alpha=0.5)
plt.ylim([-50, -10])
# plt.xticks(range(2), [])
plt.xlabel(None)
plt.ylabel('E/E\'')
plt.savefig(os.path.join(img_out, f'HFpEF_E_E_prime.{FMT}'), **fig_kwargs)


