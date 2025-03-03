import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'

ds_x = 'HEK293_ctrl_RTA'
ds_y = 'HEK293_TRUB1_kd_RTA'

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

plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
xylim = [-5, 105]
plt.scatter(df_merged_filtered['frequency_x'], df_merged_filtered['frequency_y'], s=1)
plt.plot([0, 100], [0, 100], 'r--')
plt.xlim(xylim)
plt.ylim(xylim)
plt.xlabel(ds_x, fontsize=12)
plt.ylabel(ds_y, fontsize=12)
plt.title(f'{num_sites} ${dict_mod_display[comp_mod_name]}$ sites')
# plt.savefig(os.path.join(img_out, f'scatter_{comp_mod_name}_{ds_x}_{ds_y}.png'), bbox_inches='tight')

plt.subplot(1, 2, 2)
plt.hist(df_merged_filtered['delta'], range=[-100, 100], bins=100)
plt.yscale('log')
plt.axvline(x=0, c='r', ls='--')
plt.xlabel('$\Delta$S', fontsize=12)
plt.ylabel('Site count', fontsize=12)
plt.title(f'{ds_y} - {ds_x}')

plt.savefig(os.path.join(img_out, f'{ds_x}_{ds_y}_{comp_mod_name}.png'), bbox_inches='tight')
