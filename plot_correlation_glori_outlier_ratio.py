import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

cluster_thresh = 0.8
df_file = '/home/adrian/Data/GLORI/df_outlier_ratios_thresh{:.2f}.tsv'.format(cluster_thresh)
df = pd.read_csv(df_file, sep='\t')
corr = np.corrcoef(df['Ratio'], df['ratio_outlier'])[0, 1]

plt.figure(figsize=(6, 6))
plt.plot(df['Ratio'], df['ratio_outlier'], 'o', mfc='none')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('GLORI')
plt.ylabel('Outlier')
plt.title('Cluster thresh = {:.2f}\nC = {:.2f}'.format(cluster_thresh, corr))