import os
import pandas as pd
import numpy as np

prj = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian'

mafia_file = os.path.join(prj, 'results/mAFiA_train_ISA_test_100_WT_0_IVT_sites.tsv')
m6anet_file = os.path.join(prj, 'm6Anet/100_WT_0_IVT/data.site_proba.csv.glori')
cheui_file = os.path.join(prj, 'CHEUI/100_WT_0_IVT/site_level_m6A_predictions.txt.glori')

P_VAL_THRESH = 1.0E-3
COV_THRESH = 50
def filter_df(in_df):
    df_filtered = in_df[
        (np.float32(in_df['P_adjust']) <= P_VAL_THRESH)
        & (in_df['num_features']>=COV_THRESH)
        # & (df['pred_motif']==sel_motif)
    ]
    return df_filtered

df_mafia = pd.read_csv(mafia_file, sep='\t')
df_m6anet = pd.read_csv(m6anet_file, sep=',')
df_m6anet = df_m6anet.drop(columns=['Unnamed: 0'])
df_cheui = pd.read_csv(cheui_file, sep='\t', index_col=False)

def collapse_isoforms_m6anet(in_df):
    genomic_sites = list(set([(row[0], row[1]) for row in in_df[['Chr', 'Sites']].values]))

    out_df = pd.DataFrame()
    for chr, pos in genomic_sites:
        sub_df = in_df[(in_df['Chr']==chr) & (in_df['Sites']==pos)]
        if len(sub_df)>1:
            all_mod_ratio = sub_df['mod_ratio'].values
            all_n_reads = sub_df['n_reads'].values
            total_n_reads = np.sum(all_n_reads)
            weighted_mod_ratio = np.sum(all_mod_ratio * all_n_reads) / total_n_reads
            new_row = sub_df.drop(columns=['transcript_id', 'transcript_position', 'probability_modified']).iloc[0, :]
            new_row['n_reads'] = total_n_reads
            new_row['mod_ratio'] = weighted_mod_ratio
            out_df = pd.concat([out_df, pd.DataFrame([new_row])])
        else:
            out_df = pd.concat([out_df, sub_df.drop(columns=['transcript_id', 'transcript_position', 'probability_modified'])])
    return out_df

def collapse_isoforms_cheui(in_df):
    genomic_sites = list(set([(row[0], row[1]) for row in in_df[['Chr', 'Sites']].values]))

    out_df = pd.DataFrame()
    for chr, pos in genomic_sites:
        sub_df = in_df[(in_df['Chr']==chr) & (in_df['Sites']==pos)]
        if len(sub_df)>1:
            all_mod_ratio = sub_df['stoichiometry'].values
            all_n_reads = sub_df['coverage'].values
            total_n_reads = np.sum(all_n_reads)
            weighted_mod_ratio = np.sum(all_mod_ratio * all_n_reads) / total_n_reads
            new_row = sub_df.drop(columns=['contig', 'position', 'probability', 'Strand']).iloc[0, :]
            new_row['coverage'] = total_n_reads
            new_row['stoichiometry'] = weighted_mod_ratio
            out_df = pd.concat([out_df, pd.DataFrame([new_row])])
        else:
            out_df = pd.concat([out_df, sub_df.drop(columns=['contig', 'position', 'probability', 'Strand'])])
    return out_df

df_mafia_filtered = filter_df(df_mafia)
df_m6anet_collapsed = collapse_isoforms_m6anet(df_m6anet).rename(columns={'Pvalue': 'P_adjust', 'n_reads': 'num_features'})
# df_m6anet_filtered = filter_df(df_m6anet_collapsed)
df_cheui_collapsed = collapse_isoforms_cheui(df_cheui).rename(columns={'Pvalue': 'P_adjust', 'coverage': 'num_features'})
# df_cheui_filtered = filter_df(df_cheui_collapsed)

common_cols = ['Chr', 'Sites', 'Ratio']
df_merge1 = pd.merge(df_mafia_filtered, df_m6anet_collapsed, on=common_cols, suffixes=('_mafia', '_m6anet'))
df_merge2 = pd.merge(df_merge1, df_cheui_collapsed, on=common_cols, suffixes=('', '_cheui'))
df_merge2.to_csv(os.path.join(prj, 'results/merged_comparison_mAFiA_m6Anet_CHEUI.tsv'), index=False, sep='\t')

corr_mafia = np.corrcoef(df_merge2[['Ratio', 'mod_ratio_mafia']].values.T)[0, 1]
corr_m6anet = np.corrcoef(df_merge2[['Ratio', 'mod_ratio_m6anet']].values.T)[0, 1]
corr_cheui = np.corrcoef(df_merge2[['Ratio', 'stoichiometry']].values.T)[0, 1]
print(f'{len(df_merge2)} common sites\nmAFiA: {corr_mafia:.2f}\nCHEUI: {corr_cheui:.2f}\nm6Anet: {corr_m6anet:.2f}')