import anndata
import popv
import pandas as pd
import numpy as np
import scanpy as sc

from scarches.dataset.trvae.data_handling import remove_sparsity
from popv import _check_nonnegative_integers

ref_dm_raw = sc.read(r"/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/ref_dm_RawCount.h5ad")
adata_filter = sc.read(r"/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/c9_aggr_ex_b1d8_prelim_filter.h5ad")
# adata_filter = remove_sparsity(adata_filter)

ref_dm_raw.obs["tech"] = "Reference"
adata_filter.obs["tech"] = "Query"

adata = ref_dm_raw.concatenate(adata_filter)

adata.layers["counts"] = adata.X.copy()

query_adata = adata[adata.obs["tech"] == "Query"]
# query_adata.obs_names_make_unique()
# query_adata.var_names = query_adata.var_names.str.upper()
# query_adata.var_names_make_unique()
assert _check_nonnegative_integers(query_adata.X) == True, 'Make sure query_adata.X contains raw_counts'
query_batch_key = 'batch'
# methods = ['bbknn', 'scvi', 'scanvi', 'svm', 'rf', 'onclass', 'scanorama']
methods = ['bbknn', 'scvi', 'scanvi', 'svm', 'rf', 'scanorama']
query_labels_key = None
unknown_celltype_label = 'unknown'

ref_adata = adata[adata.obs["tech"] == "Reference"]
# ref_adata.obs_names_make_unique()
# ref_adata.var_names = ref_adata.var_names.str.upper()
# ref_adata.var_names_make_unique()
assert _check_nonnegative_integers(ref_adata.X) == True, 'Make sure ref_adata.X contains raw_counts'
ref_labels_key = 'cell_type'
ref_batch_key = 'batch'

min_celltype_size = np.min(ref_adata.obs.groupby('cell_type').size())
n_samples_per_label = np.max((min_celltype_size, 100))

# import scvi

# ref_adata.layers["counts"] = ref_adata.X.copy()
# sc.pp.normalize_total(ref_adata, target_sum=1e4)
# sc.pp.log1p(ref_adata)
# ref_adata.raw = ref_adata
# sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, layer="counts", batch_key="tech", subset=True)

# scvi.model.SCVI.setup_anndata(ref_adata, layer="counts", batch_key="tech")
# vae = scvi.model.SCVI(ref_adata, n_layers=2, n_latent=30)
# vae.train()

# ref_adata.obsm["X_scVI"] = vae.get_latent_representation()

# from scvi.model.utils import mde
# ref_adata.obsm["X_mde"] = mde(ref_adata.obsm["X_scVI"])

from popv import process_query

adata = process_query(
    query_adata,
    ref_adata,
    save_folder="/home/sdxucl/Desktop/popV_ex_b1d8",
    query_batch_key=query_batch_key,
    query_labels_key=query_labels_key,
    unknown_celltype_label=unknown_celltype_label,
    pretrained_scvi_path=None,
    ref_labels_key=ref_labels_key,
    ref_batch_key=ref_batch_key,
    n_samples_per_label=n_samples_per_label,
    training_mode="offline"
)

# results_file = '/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/popv_process.h5ad'
# adata.write(results_file)

popv.annotate_data(
    adata,
    methods,
    save_path="/home/sdxucl/Desktop/popV_lessMethods",
    pretrained_scvi_path=None,
    pretrained_scanvi_path=None,
    onclass_ontology_file="/home/sdxucl/Desktop/popV/ontology/cl.ontology",
    onclass_obo_fp="/home/sdxucl/Desktop/popV/ontology/cl.obo",
    onclass_emb_fp="/home/sdxucl/Desktop/popV/ontology/cl.ontology.nlp.emb",
)

adata.write('/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/popv_ex_b1d8.h5ad')

# adata = sc.read('/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/popv_process_stuck.h5ad')

predictions = adata.obs[["popv_knn_on_scvi_offline_prediction",
        "popv_scanvi_offline_prediction",
        "popv_svm_prediction",
        "popv_rf_prediction",
        "popv_knn_on_scanorama_prediction",
        "popv_majority_vote_prediction"]]

for col in predictions.columns:
    query_adata.obs[col] = predictions.loc[query_adata.obs_names][col]

from popv import make_agreement_plots

all_prediction_keys = ["popv_knn_on_scvi_offline_prediction",
        "popv_scanvi_offline_prediction",
        "popv_svm_prediction",
        "popv_rf_prediction",
        "popv_knn_on_scanorama_prediction"]

obs_keys = query_adata.obs.keys()
pred_keys = [key for key in obs_keys if key in all_prediction_keys]

make_agreement_plots(query_adata,
                     methods=pred_keys,
                     popv_prediction_key='popv_majority_vote_prediction',
                     save_folder="/home/sdxucl/Desktop/popV_ex_b1d8")

query_adata.obs['bc2'] = query_adata.obs.index.map(lambda x: x[:-2])
cell_mapper_scvi = dict(zip(query_adata.obs.bc2, query_adata.obs.popv_knn_on_scvi_offline_prediction))
cell_mapper_scanvi = dict(zip(query_adata.obs.bc2, query_adata.obs.popv_scanvi_offline_prediction))
cell_mapper_svm = dict(zip(query_adata.obs.bc2, query_adata.obs.popv_svm_prediction))
cell_mapper_rf = dict(zip(query_adata.obs.bc2, query_adata.obs.popv_rf_prediction))
cell_mapper_scanorama = dict(zip(query_adata.obs.bc2, query_adata.obs.popv_knn_on_scanorama_prediction))
cell_mapper_popv = dict(zip(query_adata.obs.bc2, query_adata.obs.popv_majority_vote_prediction))

corr_data = sc.read('c9_aggr_ex_b1d8_mnn_correction.h5ad')
adata_filter = sc.read('c9_aggr_ex_b1d8_prelim_filter.h5ad')
sc.pp.normalize_total(adata_filter, target_sum=1e4)
sc.pp.log1p(adata_filter)
corr_data.raw = adata_filter

corr_data.obs.index = corr_data.obs.index.map(lambda x: x[:-2])
corr_data.obs["cell_type_scvi"] = corr_data.obs.index.map(cell_mapper_scvi)
corr_data.obs["cell_type_scanvi"] = corr_data.obs.index.map(cell_mapper_scanvi)
corr_data.obs["cell_type_svm"] = corr_data.obs.index.map(cell_mapper_svm)
corr_data.obs["cell_type_rf"] = corr_data.obs.index.map(cell_mapper_rf)
corr_data.obs["cell_type_scanorama"] = corr_data.obs.index.map(cell_mapper_scanorama)
corr_data.obs["cell_type_popv"] = corr_data.obs.index.map(cell_mapper_popv)

sc.pp.neighbors(corr_data, n_neighbors=20, n_pcs=40)
sc.tl.leiden(corr_data, key_added='clusters', resolution=10.0)

sc.pl.umap(corr_data, color='clusters', s=10, frameon=False, palette=color_list, legend_fontsize="small")
sc.pl.umap(corr_data, color='cell_type_scvi', s=10, frameon=False, palette=color_list, legend_fontsize="x-small",
           legend_loc="on data", legend_fontoutline=1, legend_fontweight="normal")
sc.pl.umap(corr_data, color='cell_type_scanvi', s=10, frameon=False, palette=color_list, legend_fontsize="x-small",
           legend_loc="on data", legend_fontoutline=1, legend_fontweight="normal")
sc.pl.umap(corr_data, color='cell_type_svm', s=10, frameon=False, palette=color_list, legend_fontsize="x-small",
           legend_loc="on data", legend_fontoutline=1, legend_fontweight="normal")
sc.pl.umap(corr_data, color='cell_type_rf', s=10, frameon=False, palette=color_list, legend_fontsize="x-small",
           legend_loc="on data", legend_fontoutline=1, legend_fontweight="normal")
sc.pl.umap(corr_data, color='cell_type_scanorama', s=10, frameon=False, palette=color_list, legend_fontsize="x-small",
           legend_loc="on data", legend_fontoutline=1, legend_fontweight="normal")
sc.pl.umap(corr_data, color='cell_type_popv', s=10, frameon=False, palette=color_list, legend_fontsize="x-small",
           legend_loc="on data", legend_fontoutline=1, legend_fontweight="normal")

corr_data.write("corr_data_ex_b1d8_popv_annotated.h5ad")

df_clusterOVERcelltype = corr_data.obs.groupby(["clusters", "cell_type_popv"]).size().unstack(fill_value=0)
conf_mat_clusterOVERcelltype = df_clusterOVERcelltype / df_clusterOVERcelltype.sum(axis=1).values[:, np.newaxis]
conf_mat_clusterOVERcelltype.to_csv("/home/sdxucl/Desktop/popV_ex_b1d8/conf_mat_clusterOVERcelltype.csv")
df_celltypeOVERcluster = corr_data.obs.groupby(["cell_type_popv", "clusters"]).size().unstack(fill_value=0)
conf_celltypeOVERcluster = df_celltypeOVERcluster / df_celltypeOVERcluster.sum(axis=1).values[:, np.newaxis]
conf_celltypeOVERcluster.to_csv("/home/sdxucl/Desktop/popV_ex_b1d8/conf_mat_celltypeOVERcluster.csv")

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 8))
_ = plt.pcolor(conf_mat)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90, size=5)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index, size=5)
plt.xlabel("cell_type_popv", size=15)
plt.ylabel("clusters", size=15)
plt.savefig("/home/sdxucl/Desktop/scvi_label_transfer/con_mx.png", dpi=800)
