import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import mnnpy

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

sc.settings.set_figure_params(dpi=200, facecolor='white')

results_file_MnnCorrection = '/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/' \
                             'c9_aggr_mnn_correction.h5ad'
results_file_prelim_filter = '/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/' \
                             'c9_aggr_prelim_filter.h5ad'

print("Reading data...")
adata = sc.read_10x_mtx(
    '/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/outs/count/filtered_feature_bc_matrix',
    var_names='gene_symbols', cache=True)

adata.var_names_make_unique()

adata.obs['p_number'] = adata.obs_names.str.replace(r'.*-', '').astype('category')

adata.obs.loc[adata.obs['p_number'].isin(['1', '2', '3', '4', '5']), 'batch'] = '1'
adata.obs.loc[adata.obs['p_number'].isin(['6', '7', '8', '9', '10']), 'batch'] = '2'
adata.obs.loc[adata.obs['p_number'].isin(['11', '12', '13', '14', '15']), 'batch'] = '3'

adata.obs.loc[adata.obs['p_number'].isin(['1', '6', '11']), 'day'] = '0'
adata.obs.loc[adata.obs['p_number'].isin(['2', '7', '12']), 'day'] = '2'
adata.obs.loc[adata.obs['p_number'].isin(['3', '8', '13']), 'day'] = '4'
adata.obs.loc[adata.obs['p_number'].isin(['4', '9', '14']), 'day'] = '6'
adata.obs.loc[adata.obs['p_number'].isin(['5', '10', '15']), 'day'] = '8'

adata.obs.loc[adata.obs['p_number'] == '1', 'sample'] = 'batch1_day0'
adata.obs.loc[adata.obs['p_number'] == '2', 'sample'] = 'batch1_day2'
adata.obs.loc[adata.obs['p_number'] == '3', 'sample'] = 'batch1_day4'
adata.obs.loc[adata.obs['p_number'] == '4', 'sample'] = 'batch1_day6'
adata.obs.loc[adata.obs['p_number'] == '5', 'sample'] = 'batch1_day8'
adata.obs.loc[adata.obs['p_number'] == '6', 'sample'] = 'batch2_day0'
adata.obs.loc[adata.obs['p_number'] == '7', 'sample'] = 'batch2_day2'
adata.obs.loc[adata.obs['p_number'] == '8', 'sample'] = 'batch2_day4'
adata.obs.loc[adata.obs['p_number'] == '9', 'sample'] = 'batch2_day6'
adata.obs.loc[adata.obs['p_number'] == '10', 'sample'] = 'batch2_day8'
adata.obs.loc[adata.obs['p_number'] == '11', 'sample'] = 'batch3_day0'
adata.obs.loc[adata.obs['p_number'] == '12', 'sample'] = 'batch3_day2'
adata.obs.loc[adata.obs['p_number'] == '13', 'sample'] = 'batch3_day4'
adata.obs.loc[adata.obs['p_number'] == '14', 'sample'] = 'batch3_day6'
adata.obs.loc[adata.obs['p_number'] == '15', 'sample'] = 'batch3_day8'

adata = adata[adata.obs['p_number'].isin(['1', '2', '3', '4', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'])]  # take out "5" as an outlier

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('mt:')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

adata_filtered = adata
# adata_filtered.write_h5ad(results_file_prelim_filter)


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

print(adata.X.shape)

adata2 = adata.raw.to_adata()
print(adata2.X[1:10, 1:10])

var_genes_all = adata.var.highly_variable
print("Highly variable genes: %d" % sum(var_genes_all))

sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='batch')
print("Highly variable genes intersection: %d" % sum(adata2.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata2.var.highly_variable_nbatches.value_counts())

var_genes_batch = adata2.var.highly_variable_nbatches > 0

print("Any batch var genes: %d" % sum(var_genes_batch))
print("All data var genes: %d" % sum(var_genes_all))
print("Overlap: %d" % sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d" % sum(adata2.var.highly_variable_nbatches == 3))
print("Overlap batch instersection and all: %d" % sum(var_genes_all & adata2.var.highly_variable_intersection))

var_select = adata2.var.highly_variable_nbatches > 2
var_genes = var_select.index[var_select]
len(var_genes)

adata.obs['batch'] = pd.Categorical(adata.obs.batch)
batches = adata.obs['batch'].cat.categories.tolist()
alldata = {}
for batch in batches:
    alldata[batch] = adata2[adata2.obs['batch'] == batch,]

alldata

import numba as nb

cdata = mnnpy.mnn_correct(alldata['1'], alldata['2'], alldata['3'], var_subset=var_genes)
corr_data = cdata[0][:, var_genes]
corr_data.X.shape

sc.tl.pca(corr_data, svd_solver='arpack', use_highly_variable=False)
sc.pl.pca(corr_data, components=['1,2', '3,4', '5,6', '7,8'], ncols=2, color='batch')

sc.pp.neighbors(adata, n_pcs=40, n_neighbors=20)
sc.tl.umap(adata)

sc.pp.neighbors(corr_data, n_pcs=40, n_neighbors=20)
sc.tl.umap(corr_data)

sc.pl.umap(corr_data, color='batch', title='MNN Corrected UMAP', frameon=False, s=5)
sc.pl.umap(corr_data, color='day', title='MNN Corrected UMAP', frameon=False, s=5)
sc.pl.umap(adata, color='batch', title='Uncorrected UMAP', frameon=False)
sc.pl.umap(adata, color='day', title='Uncorrected UMAP', frameon=False, s=5)
sc.pl.umap(adata, color='sample', title='Uncorrected UMAP', frameon=False, s=5)

cdata1 = cdata[0][:, var_genes]

corr_data_batch1 = corr_data[corr_data.obs['batch'].isin(['0'])]
corr_data_batch2 = corr_data[corr_data.obs['batch'].isin(['1'])]
corr_data_batch3 = corr_data[corr_data.obs['batch'].isin(['2'])]

sc.pl.umap(corr_data_batch1, color='day', s=5, frameon=False, title='Batch1')
sc.pl.umap(corr_data_batch2, color='day', s=5, frameon=False, title='Batch2')
sc.pl.umap(corr_data_batch3, color='day', s=5, frameon=False, title='Batch3')

adata_filtered.write_h5ad('/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/c9_aggr_ex_b1d8_prelim_filter.h5ad')
corr_data.write_h5ad('/home/sdxucl/Work/CountedMatrices/C9_AGGR_ALL_NO_Normalization/c9_aggr_ex_b1d8_mnn_correction.h5ad')
