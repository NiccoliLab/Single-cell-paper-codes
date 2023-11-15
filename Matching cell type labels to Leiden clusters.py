import numpy as np
import pandas as pd

conf_mat_clusterOVERcelltype = pd.read_csv("/home/sdxucl/Desktop/Old_Desktop/popV/conf_mat_clusterOVERcelltype.csv")
# conf_celltypeOVERcluster = pd.read_csv("/home/sdxucl/Desktop/popV_ex_b1d8/conf_mat_celltypeOVERcluster.csv")

df = conf_mat_clusterOVERcelltype.T

df = df.drop(['clusters'])

a = [0, df[0].idxmax(), df[0].max()]

df_cluster_popv_label = pd.DataFrame(a)

for i in range(1, 214):
    df_cluster_popv_label[i] = [i, df[i].idxmax(), df[i].max()]

df_cluster_popv_label = df_cluster_popv_label.T
df_cluster_popv_label = df_cluster_popv_label.rename(columns={0: "leiden_clusters", 1: "popV_label", 2: "percentage"})
df_cluster_popv_label = df_cluster_popv_label.set_index("popV_label")
df_cluster_popv_label = df_cluster_popv_label.sort_index()
# df_cluster_popv_label = df_cluster_popv_label.sort_values(by=['percentage'], ascending=False)
df_cluster_popv_label.to_csv("/home/sdxucl/Desktop/popV/df_cluster_popv_label.csv")

leiden_label_over90 = df_cluster_popv_label[df_cluster_popv_label['percentage'] >= 0.9]
leiden_label_over90['popV_label'] = leiden_label_over90.index
leiden_label_over90 =leiden_label_over90.set_index('leiden_clusters')
del leiden_label_over90['percentage']
g = leiden_label_over90.groupby('popV_label')
leiden_label_over90.loc[g['popV_label'].transform('size').gt(1), 'popV_label'] += '_subtype' + g.cumcount().astype(str)
leiden_label_over90['popV_label'] = 'popV_' + leiden_label_over90['popV_label'].astype(str)
leiden_label_over90.to_csv("/home/sdxucl/Desktop/popV/leiden_label_over90.csv")

leiden_label_over80 = df_cluster_popv_label[df_cluster_popv_label['percentage'] >= 0.8]
leiden_label_over80['popV_label'] = leiden_label_over80.index
leiden_label_over80 =leiden_label_over80.set_index('leiden_clusters')
del leiden_label_over80['percentage']
g = leiden_label_over80.groupby('popV_label')
leiden_label_over80.loc[g['popV_label'].transform('size').gt(1), 'popV_label'] += '_subtype' + g.cumcount().astype(str)
leiden_label_over80['popV_label'] = 'popV_' + leiden_label_over80['popV_label'].astype(str)
leiden_label_over80.to_csv("/home/sdxucl/Desktop/popV/leiden_label_over80.csv")

leiden_label_over70 = df_cluster_popv_label[df_cluster_popv_label['percentage'] >= 0.7]
leiden_label_over70['popV_label'] = leiden_label_over70.index
leiden_label_over70 =leiden_label_over70.set_index('leiden_clusters')
del leiden_label_over70['percentage']
g = leiden_label_over70.groupby('popV_label')
leiden_label_over70.loc[g['popV_label'].transform('size').gt(1), 'popV_label'] += '_subtype' + g.cumcount().astype(str)
leiden_label_over70['popV_label'] = 'popV_' + leiden_label_over70['popV_label'].astype(str)
leiden_label_over70.index = leiden_label_over70.index.astype(str)
leiden_label_over70.to_csv("/home/sdxucl/Desktop/popV/leiden_label_over70.csv")

dict_popv_clusters = leiden_label_over70.to_dict('dict')
dict_popv_clusters['popV_label']

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['figure.figsize'] = (8, 8)  # rescale figures
sc.settings.verbosity = 3

import milopy
import milopy.core as milo
import milopy.plot as milopl

corr_data = sc.read('c9_aggr_mnn_correction.h5ad')
adata_filter = sc.read('c9_aggr_prelim_filter.h5ad')
sc.pp.normalize_total(adata_filter, target_sum=1e4)
sc.pp.log1p(adata_filter)
corr_data.raw = adata_filter

sc.pp.neighbors(corr_data, n_neighbors=20, n_pcs=40)
sc.tl.leiden(corr_data, key_added='clusters', resolution=10.0)

corr_data.obs = corr_data.obs.replace({"clusters": dict_popv_clusters['popV_label']})

sc.pl.umap(corr_data, color='clusters', s=50, frameon=False, palette=color_list, legend_fontsize="small", legend_loc='on data', legend_fontoutline=1)
