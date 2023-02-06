import os
HOME = os.path.expanduser('~')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

cluster_thresh = 0.8
df_file = '/home/adrian/Data/GLORI/df_outlier_ratios_thresh{:.2f}.tsv'.format(cluster_thresh)
df_all = pd.read_csv(df_file, sep='\t')
# df = df_all[df_all['P_adjust']==0]
df = df_all
corr = np.corrcoef(df['Ratio'], df['ratio_outlier'])[0, 1]

img_out = os.path.join(HOME, 'img_out')

plt.figure(figsize=(6, 6))
plt.plot(df['NormeRatio'], df['ratio_outlier'], 'o', mfc='none')
plt.plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1), 'r-', alpha=0.5)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('GLORI')
plt.ylabel('Outlier')
plt.title('Cluster thresh = {:.2f}\nC = {:.2f}'.format(cluster_thresh, corr))
plt.show()

# plt.savefig(os.path.join(img_out, 'corr_glori_outlier.png'), bbox_inches='tight')
# plt.close()