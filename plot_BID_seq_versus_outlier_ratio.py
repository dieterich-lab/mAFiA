import pandas as pd
import matplotlib.pyplot as plt

mod_file = '/home/adrian/img_out/MAFIA_clustering/BID_seq/res_outlier_ratio_MaxAbs.tsv'
df_mod = pd.read_csv(mod_file, sep='\t')

plt.figure(figsize=(5, 5))
plt.scatter(df_mod['Frac_Ave %'], df_mod['outlier_ratio'], s=50, facecolors='none', edgecolors='blue')
plt.xlim([-5, 105])
plt.ylim([-5, 105])
plt.xlabel('BID-Seq Fraction Avg.', fontsize=12)
plt.ylabel('ONT Outlier Fraction', fontsize=12)
plt.title('$\psi$ mod in HEK293A WT\n{} Sites with $\geq$50 reads'.format(len(df_mod)), fontsize=15)

plt.savefig('/home/adrian/img_out/MAFIA_clustering/BID_seq/BID_seq_frac_versus_outlier_ratio_MaxAbs.png', bbox_inches='tight')
plt.close('all')