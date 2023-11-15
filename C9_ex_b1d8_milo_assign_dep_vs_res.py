# %load_ext autoreload
# %autoreload 2

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

corr_data = sc.read('c9_aggr_ex_b1d8_mnn_correction.h5ad')
adata_filter = sc.read('c9_aggr_ex_b1d8_prelim_filter.h5ad')
sc.pp.normalize_total(adata_filter, target_sum=1e4)
sc.pp.log1p(adata_filter)
corr_data.raw = adata_filter

sc.pp.neighbors(corr_data, n_neighbors=20, n_pcs=40)
sc.tl.leiden(corr_data, key_added='clusters', resolution=10.0)
sc.pl.umap(corr_data, color='clusters', s=50, frameon=False, palette=color_list, legend_fontsize="small",
           legend_loc='on data', legend_fontoutline=1)

# corr_data_0day = corr_data[corr_data.obs['day'] == '0']

# sc.tl.rank_genes_groups(corr_data, 'clusters', method='wilcoxon')
# sc.pl.rank_genes_groups(corr_data_0day, n_genes=25)

corr_data.obs.groupby("clusters").apply(len)

# marker_genes = sc.get.rank_genes_groups_df(corr_data, group=None)
# marker_genes.to_csv("/home/sdxucl/Desktop/MarkerGenesCluster10.0.csv")

milo.make_nhoods(corr_data, prop=0.1)

corr_data.obsm["nhoods"]

corr_data[corr_data.obs['nhood_ixs_refined'] != 0].obs[['nhood_ixs_refined', 'nhood_kth_distance']]

nhood_size = np.array(corr_data.obsm["nhoods"].sum(0)).ravel()
plt.hist(nhood_size, bins=100);  # adjust the value of n_neighbors to achieve number of cells per neighbourhood to
# peak at least 3 x number of groups cells

milo.count_nhoods(corr_data, sample_col="sample")

corr_data.uns["nhood_adata"]

corr_data.obs["days_continuous"] = corr_data.obs["day"].cat.codes

milo.DA_nhoods(corr_data, design="~ days_continuous")

corr_data.uns["nhood_adata"].obs  # to_csv(r"/home/sdxucl/Desktop/check.csv")

old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [10, 5]
plt.subplot(1, 2, 1)
plt.hist(corr_data.uns["nhood_adata"].obs.PValue, bins=50);
plt.xlabel("P-Vals");
plt.subplot(1, 2, 2)
plt.plot(corr_data.uns["nhood_adata"].obs.logFC, -np.log10(corr_data.uns["nhood_adata"].obs.SpatialFDR), '.');
plt.xlabel("log-Fold Change");
plt.ylabel("- log10(Spatial FDR)");
plt.tight_layout()
plt.rcParams["figure.figsize"] = old_figsize

import milopy.utils

milopy.utils.build_nhood_graph(corr_data)

plt.rcParams["figure.figsize"] = [10, 10]
plt.rcParams.update({'font.size': 22})
ax = milopl.plot_nhood_graph(corr_data,
                        alpha=0.1,  ## SpatialFDR level (10%)
                        min_size=5,  ## Size of the smallest dot
                        )
plt.suptitle(f'Milo Differential Abundance Analysis', size=50)

# corr_data.uns["nhood_adata"].obs["logFC"].values[corr_data.uns["nhood_adata"].obs["logFC"].values > 2] = 2
# corr_data.uns["nhood_adata"].obs["logFC"].values[corr_data.uns["nhood_adata"].obs["logFC"].values < -2] = -2

# milopy.utils.annotate_nhoods(corr_data, anno_col='clusters')
milopy.utils.annotate_nhoods(corr_data, anno_col='clusters')

plt.hist(corr_data.uns['nhood_adata'].obs["nhood_annotation_frac"]);
plt.xlabel("clusters fraction")

sc.pl.violin(corr_data.uns['nhood_adata'], "logFC", groupby="nhood_annotation", rotation=90, show=False,
             palette="prism");
plt.axhline(y=0, color='black', linestyle='--');
plt.show()

df_nhood_FDR10 = corr_data.uns['nhood_adata'].obs[corr_data.uns['nhood_adata'].obs['FDR'] < 0.1]
df_nhood_FDR10['changed_cell_number'] = (2 ** df_nhood_FDR10['logFC'] - 1) * df_nhood_FDR10['Nhood_size'] * \
                                        df_nhood_FDR10['nhood_annotation_frac']

df_clusters_nhood_FDR10 = df_nhood_FDR10.groupby(df_nhood_FDR10["nhood_annotation"]).sum()
df_clusters_nhood_FDR10 = df_clusters_nhood_FDR10['changed_cell_number']
df_clusters_nhood_FDR10 = df_clusters_nhood_FDR10.to_frame()
df_clusters_nhood_FDR10['nhood_annotation'] = df_clusters_nhood_FDR10.index

cell_num_clusters = corr_data.obs.groupby("clusters").apply(len)
cell_num_clusters = cell_num_clusters.to_frame(name="cell_numbers")
cell_num_clusters["nhood_annotation"] = cell_num_clusters.index

df_clusters_nhood_FDR10.reset_index(drop=True, inplace=True)
cell_num_clusters.reset_index(drop=True, inplace=True)
df_clusters_nhood_FDR10 = df_clusters_nhood_FDR10.merge(cell_num_clusters, how="inner", on="nhood_annotation")
df_clusters_nhood_FDR10["normalized_changed_cell_num"] = df_clusters_nhood_FDR10["changed_cell_number"] / \
                                                         df_clusters_nhood_FDR10["cell_numbers"]
df_clusters_nhood_FDR10 = df_clusters_nhood_FDR10.loc[~(df_clusters_nhood_FDR10["normalized_changed_cell_num"] == 0)]

ax = df_clusters_nhood_FDR10.plot.scatter(x="nhood_annotation", y="normalized_changed_cell_num",
                                          c="normalized_changed_cell_num", colormap='viridis')
plt.axhline(y=0, color='black', linestyle='--');
plt.axhline(y=0.1, color='black', linestyle='--');
plt.axhline(y=-0.1, color='black', linestyle='--');
plt.show()

nhood_FDR10 = corr_data.uns['nhood_adata'][corr_data.uns['nhood_adata'].obs['FDR'] <= 0.1]
sc.pl.violin(nhood_FDR10, "logFC", groupby="nhood_annotation", rotation=90, show=False, palette="prism");
plt.axhline(y=0, color='black', linestyle='--');
plt.show()

# df_nhood_FDR15 = corr_data.uns['nhood_adata'].obs[corr_data.uns['nhood_adata'].obs['FDR'] < 0.15]
# df_nhood_FDR15['changed_cell_number'] = (2 ** df_nhood_FDR15['logFC'] - 1) * df_nhood_FDR15['Nhood_size'] * df_nhood_FDR15['nhood_annotation_frac']

# df_clusters_nhood_FDR15 = df_nhood_FDR15.groupby(df_nhood_FDR15["nhood_annotation"]).sum()
# df_clusters_nhood_FDR15 = df_clusters_nhood_FDR15['changed_cell_number']
# df_clusters_nhood_FDR15 = df_clusters_nhood_FDR15.to_frame()
# df_clusters_nhood_FDR15['nhood_annotation'] = df_clusters_nhood_FDR15.index

# df_clusters_nhood_FDR15.reset_index(drop=True, inplace=True)
# df_clusters_nhood_FDR15 = df_clusters_nhood_FDR15.merge(cell_num_clusters, how="inner", on="nhood_annotation")
# df_clusters_nhood_FDR15["normalized_changed_cell_num"] = df_clusters_nhood_FDR15["changed_cell_number"] / df_clusters_nhood_FDR15["cell_numbers"]
# df_clusters_nhood_FDR15 = df_clusters_nhood_FDR15.loc[~(df_clusters_nhood_FDR15["normalized_changed_cell_num"] == 0)]

# ax =df_clusters_nhood_FDR15.plot.scatter(x="nhood_annotation", y="normalized_changed_cell_num", c="normalized_changed_cell_num", colormap='viridis')
# plt.axhline(y=0, color='black', linestyle='--');
# plt.axhline(y=0.1, color='black', linestyle='--');
# plt.axhline(y=-0.1, color='black', linestyle='--');
# plt.show()

# nhood_FDR15 = corr_data.uns['nhood_adata'][corr_data.uns['nhood_adata'].obs['FDR'] < 0.15]
# sc.pl.violin(nhood_FDR15, "logFC", groupby="nhood_annotation", rotation=90, show=False, palette="prism");
# plt.axhline(y=0, color='black', linestyle='--');
# plt.show()

# df_nhood_FDR20 = corr_data.uns['nhood_adata'].obs[corr_data.uns['nhood_adata'].obs['FDR'] < 0.2]
# df_nhood_FDR20['changed_cell_number'] = (2 ** df_nhood_FDR20['logFC'] - 1) * df_nhood_FDR20['Nhood_size'] * df_nhood_FDR20['nhood_annotation_frac']

# df_clusters_nhood_FDR20 = df_nhood_FDR20.groupby(df_nhood_FDR20["nhood_annotation"]).sum()
# df_clusters_nhood_FDR20 = df_clusters_nhood_FDR20['changed_cell_number']
# df_clusters_nhood_FDR20 = df_clusters_nhood_FDR20.to_frame()
# df_clusters_nhood_FDR20['nhood_annotation'] = df_clusters_nhood_FDR20.index

# df_clusters_nhood_FDR20.reset_index(drop=True, inplace=True)
# df_clusters_nhood_FDR20 = df_clusters_nhood_FDR20.merge(cell_num_clusters, how="inner", on="nhood_annotation")
# df_clusters_nhood_FDR20["normalized_changed_cell_num"] = df_clusters_nhood_FDR20["changed_cell_number"] / df_clusters_nhood_FDR20["cell_numbers"]
# df_clusters_nhood_FDR20 = df_clusters_nhood_FDR20.loc[~(df_clusters_nhood_FDR20["normalized_changed_cell_num"] == 0)]

# ax =df_clusters_nhood_FDR20.plot.scatter(x="nhood_annotation", y="normalized_changed_cell_num", c="normalized_changed_cell_num", colormap='viridis')
# plt.axhline(y=0, color='black', linestyle='--');
# plt.axhline(y=0.1, color='black', linestyle='--');
# plt.axhline(y=-0.1, color='black', linestyle='--');
# plt.show()

# nhood_FDR20 = corr_data.uns['nhood_adata'][corr_data.uns['nhood_adata'].obs['FDR'] < 0.2]
# sc.pl.violin(nhood_FDR20, "logFC", groupby="nhood_annotation", rotation=90, show=False, palette="prism");
# plt.axhline(y=0, color='black', linestyle='--');
# plt.show()

# nhood = corr_data.uns['nhood_adata'].obs
# nhood['changed_cell_number'] = (2 ** nhood['logFC'] - 1) * nhood['Nhood_size'] * nhood['nhood_annotation_frac']

# df_clusters_nhood = nhood.groupby(nhood["nhood_annotation"]).sum()
# df_clusters_nhood = df_clusters_nhood['changed_cell_number']
# df_clusters_nhood = df_clusters_nhood.to_frame()
# df_clusters_nhood['nhood_annotation'] = df_clusters_nhood.index

# df_clusters_nhood.reset_index(drop=True, inplace=True)
# df_clusters_nhood = df_clusters_nhood.merge(cell_num_clusters, how="inner", on="nhood_annotation")
# df_clusters_nhood["normalized_changed_cell_num"] = df_clusters_nhood["changed_cell_number"] / df_clusters_nhood["cell_numbers"]
# df_clusters_nhood = df_clusters_nhood.loc[~(df_clusters_nhood["normalized_changed_cell_num"] == 0)]

# ax =df_clusters_nhood.plot.scatter(x="nhood_annotation", y="normalized_changed_cell_num", c="normalized_changed_cell_num", colormap='viridis')
# plt.axhline(y=0, color='black', linestyle='--');
# plt.axhline(y=0.1, color='black', linestyle='--');
# plt.axhline(y=-0.1, color='black', linestyle='--');
# plt.show()

dep_clusters = df_clusters_nhood_FDR10[df_clusters_nhood_FDR10["normalized_changed_cell_num"] < -0.2]
dep_clusters = dep_clusters["nhood_annotation"].to_list()
dep_clusters = list(map(int, dep_clusters))

inc_clusters = df_clusters_nhood_FDR10[df_clusters_nhood_FDR10["normalized_changed_cell_num"] > 0.2]
inc_clusters = inc_clusters["nhood_annotation"].to_list()
inc_clusters = list(map(int, inc_clusters))

all_clusters = list(range(0, 213))

res_clusters = list(set(all_clusters) - set(dep_clusters) - set(inc_clusters))

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.cm as cm

nhood_dep = pd.DataFrame(columns=corr_data.uns['nhood_adata'].obs.columns)

for x in dep_clusters:
    nhood_x = corr_data.uns['nhood_adata'].obs[corr_data.uns['nhood_adata'].obs["nhood_annotation"] == str(x)]
    nhood_dep = pd.concat([nhood_dep, nhood_x])
    fig, ax = plt.subplots()
    ax = plt.gca()
    ax = sns.swarmplot(x=nhood_x["logFC"], ax=ax, size=12)
    ax.set(xlim=(-5, 5))
    ax.set_ylabel(str(x), rotation=0, fontsize=20, labelpad=20)
    for a in ax.get_children():
        if type(a) is matplotlib.collections.PathCollection:
            offsets = a.get_offsets()
            cmap = sns.color_palette("coolwarm", as_cmap=True)
            norm = plt.Normalize(vmin=-3, vmax=3)
            facecolors = [cmap(norm(x)) for x, y in offsets]
            a.set_color(facecolors)
            break
    fig.set_size_inches(10, 8)
    fig.tight_layout()
    plt.savefig("/home/sdxucl/Desktop/swarmplot_dep_nhoodFDR10/" + str(x) + "_milo_swarm.png")
    plt.clf()

nhood_res = pd.DataFrame(columns=corr_data.uns['nhood_adata'].obs.columns)

for x in res_clusters:
    nhood_x = corr_data.uns['nhood_adata'].obs[corr_data.uns['nhood_adata'].obs["nhood_annotation"] == str(x)]
    nhood_res = pd.concat([nhood_res, nhood_x])
    fig, ax = plt.subplots()
    ax = plt.gca()
    ax = sns.swarmplot(x=nhood_x["logFC"], ax=ax, size=12)
    ax.set(xlim=(-5, 5))
    ax.set_ylabel(str(x), rotation=0, fontsize=20, labelpad=20)
    for a in ax.get_children():
        if type(a) is matplotlib.collections.PathCollection:
            offsets = a.get_offsets()
            cmap = sns.color_palette("coolwarm", as_cmap=True)
            norm = plt.Normalize(vmin=-3, vmax=3)
            facecolors = [cmap(norm(x)) for x, y in offsets]
            a.set_color(facecolors)
            break
    fig.set_size_inches(10, 8)
    fig.tight_layout()
    plt.savefig("/home/sdxucl/Desktop/swarmplot_res_nhoodFDR10/" + str(x) + "_milo_swarm.png")
    plt.clf()

nhood_res_decreasing = nhood_res[nhood_res.changed_cell_number < 0]
df_clusters_nhood_FDR10_res_decreasing = nhood_res_decreasing.groupby(nhood_res_decreasing["nhood_annotation"]).sum()
df_clusters_nhood_FDR10_res_decreasing = df_clusters_nhood_FDR10_res_decreasing['changed_cell_number']
df_clusters_nhood_FDR10_res_decreasing = df_clusters_nhood_FDR10_res_decreasing.to_frame()
df_clusters_nhood_FDR10_res_decreasing['nhood_annotation'] = df_clusters_nhood_FDR10_res_decreasing.index
df_clusters_nhood_FDR10_res_decreasing.reset_index(drop=True, inplace=True)
df_clusters_nhood_FDR10_res_decreasing = df_clusters_nhood_FDR10_res_decreasing.merge(cell_num_clusters, how="inner",
                                                                                      on="nhood_annotation")
df_clusters_nhood_FDR10_res_decreasing["changed_ratio"] = df_clusters_nhood_FDR10_res_decreasing[
                                                              "changed_cell_number"] / \
                                                          df_clusters_nhood_FDR10_res_decreasing["cell_numbers"]
mov_clusters_dec = df_clusters_nhood_FDR10_res_decreasing[df_clusters_nhood_FDR10_res_decreasing.changed_ratio < -0.2][
    "nhood_annotation"].to_list()
mov_clusters_dec = list(map(int, mov_clusters_dec))

nhood_res_increasing = nhood_res[nhood_res.changed_cell_number > 0]
df_clusters_nhood_FDR10_res_increasing = nhood_res_increasing.groupby(nhood_res_decreasing["nhood_annotation"]).sum()
df_clusters_nhood_FDR10_res_increasing = df_clusters_nhood_FDR10_res_increasing['changed_cell_number']
df_clusters_nhood_FDR10_res_increasing = df_clusters_nhood_FDR10_res_increasing.to_frame()
df_clusters_nhood_FDR10_res_increasing['nhood_annotation'] = df_clusters_nhood_FDR10_res_increasing.index
df_clusters_nhood_FDR10_res_increasing.reset_index(drop=True, inplace=True)
df_clusters_nhood_FDR10_res_increasing = df_clusters_nhood_FDR10_res_increasing.merge(cell_num_clusters, how="inner",
                                                                                      on="nhood_annotation")
df_clusters_nhood_FDR10_res_increasing["changed_ratio"] = df_clusters_nhood_FDR10_res_increasing[
                                                              "changed_cell_number"] / \
                                                          df_clusters_nhood_FDR10_res_increasing["cell_numbers"]
mov_clusters_inc = df_clusters_nhood_FDR10_res_increasing[df_clusters_nhood_FDR10_res_increasing.changed_ratio > 0.2][
    "nhood_annotation"].to_list()
mov_clusters_inc = list(map(int, mov_clusters_dec))

mov_clusters = mov_clusters_dec + list(set(mov_clusters_inc) - set(mov_clusters_dec))

res_clusters = list(set(res_clusters) - set(mov_clusters))

pd.DataFrame(dep_clusters, columns={"clusters"}).to_csv("/home/sdxucl/Desktop/dep_clusters.csv")
pd.DataFrame(res_clusters, columns={"clusters"}).to_csv("/home/sdxucl/Desktop/res_clusters.csv")
pd.DataFrame(mov_clusters, columns={"clusters"}).to_csv("/home/sdxucl/Desktop/mov_clusters.csv")

fig, ax = plt.subplots()

vcenter = 0
vmin, vmax = -3, 3
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
colormap = cm.coolwarm

sns.scatterplot(
    x="logFC",
    y="nhood_annotation",
    data=nhood_dep,
    size=nhood_dep["Nhood_size"],
    c=nhood_dep["logFC"],
    norm=normalize,
    cmap=colormap,
    ax=ax,
    sizes=(50, 500)
)

pts = ax.collections[0]
pts.set_offsets(pts.get_offsets() + np.c_[np.zeros(len(nhood_dep)), np.random.uniform(-.2, .2, len(nhood_dep))])
ax.set(xlim=(-5, 5))
ax.axvline(vcenter, color="grey", ls="--")

nhood_dep["logFC"] = nhood_dep["logFC"].astype(float)

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(nhood_dep["logFC"])
fig.colorbar(scalarmappaple)

fig, ax = plt.subplots()

vcenter = 0
vmin, vmax = -3, 3
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
colormap = cm.coolwarm

sns.scatterplot(
    x="logFC",
    y="nhood_annotation",
    data=nhood_res,
    size=nhood_res["Nhood_size"],
    c=nhood_res["logFC"],
    norm=normalize,
    cmap=colormap,
    ax=ax,
    sizes=(50, 500)
)

pts = ax.collections[0]
pts.set_offsets(pts.get_offsets() + np.c_[np.zeros(len(nhood_res)), np.random.uniform(-.2, .2, len(nhood_res))])
ax.set(xlim=(-5, 5))
ax.axvline(vcenter, color="grey", ls="--")

nhood_res["logFC"] = nhood_res["logFC"].astype(float)

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(nhood_res["logFC"])
fig.colorbar(scalarmappaple)

nhood_res = pd.DataFrame(columns=df_nhood_filtered.columns)

for x in res_clusters:
    nhood_x = df_nhood[df_nhood["nhood_annotation"] == str(x)]
    nhood_res = pd.concat([nhood_res, nhood_x])

fig, ax = plt.subplots()

vcenter = 0
vmin, vmax = -3, 3
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
colormap = cm.coolwarm

sns.scatterplot(
    x="logFC",
    y="nhood_annotation",
    data=nhood_res,
    size=nhood_res["Nhood_size"],
    c=nhood_res["logFC"],
    norm=normalize,
    cmap=colormap,
    ax=ax,
)

pts = ax.collections[0]
pts.set_offsets(pts.get_offsets() + np.c_[np.zeros(len(nhood_res)), np.random.uniform(-.2, .2, len(nhood_res))])
ax.set(xlim=(-5, 5))
ax.axvline(vcenter, color="grey", ls="--")

nhood_res["logFC"] = nhood_res["logFC"].astype(float)

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(nhood_res["logFC"])
fig.colorbar(scalarmappaple)

dep_clusters = pd.read_csv("/home/sdxucl/Desktop/dep_clusters.csv")
dep_clusters = dep_clusters[["clusters"]].astype(str)
dep_clusters["type"] = "depleted_cluster"

res_clusters = pd.read_csv("/home/sdxucl/Desktop/res_clusters.csv")
res_clusters = res_clusters[["clusters"]].astype(str)
res_clusters["type"] = "resistant_cluster"

mov_clusters = pd.read_csv("/home/sdxucl/Desktop/mov_clusters.csv")
mov_clusters = mov_clusters[["clusters"]].astype(str)
mov_clusters["type"] = "moving_cluster"

inc_clusters = pd.read_csv("/home/sdxucl/Desktop/inc_clusters.csv")
inc_clusters = inc_clusters[["clusters"]].astype(str)
inc_clusters["type"] = "increasing cluster"

cluster_no_sig_nhood = {str(x) for x in set(range(0, 214))} - set(df_clusters_nhood_FDR10["nhood_annotation"])
cluster_no_sig_nhood = pd.DataFrame(cluster_no_sig_nhood)
cluster_no_sig_nhood = cluster_no_sig_nhood.rename(columns={0: "nhood_annotation"})
cluster_no_sig_nhood["normalized_changed_cell_num"] = 0

df_clusters_nhood_FDR10 = df_clusters_nhood_FDR10[["nhood_annotation", "normalized_changed_cell_num"]]
df_clusters_nhood = pd.concat([df_clusters_nhood_FDR10, cluster_no_sig_nhood])

cluster_label = pd.concat([dep_clusters, res_clusters, mov_clusters, inc_clusters])
cluster_label = cluster_label.rename(columns={"clusters": "nhood_annotation"})
df_clusters_nhood = df_clusters_nhood.merge(cluster_label, how="left", on="nhood_annotation")

sns.displot(df_clusters_nhood, x="normalized_changed_cell_num", hue="type", height=25, aspect=25/25, stat="probability")
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)

plt.plot(x, p, 'k', linewidth=2)
plt.show()

mean = df_clusters_nhood["normalized_changed_cell_num"].mean()
std = df_clusters_nhood["normalized_changed_cell_num"].std()

# Set up the figure and axes for the plot
fig, ax1 = plt.subplots()

# Create the normal distribution curve
x = np.linspace(mean - 4*std, mean + 4*std, 100)
y = 1/(std*np.sqrt(2*np.pi)) * np.exp(-(x-mean)**2/(2*std**2))
ax1.plot(x, y, color='black', linewidth=2)

# Create the two vertical lines for the 95% confidence interval
z = 1.96  # 95% confidence interval corresponds to z-score of 1.96
ci_min = mean - z*std  # /np.sqrt(len(df_clusters_nhood))
ci_max = mean + z*std  # /np.sqrt(len(df_clusters_nhood))
# ax1.axvline(ci_min, color='gray', linestyle='--', linewidth=2)
# ax1.axvline(ci_max, color='gray', linestyle='--', linewidth=2)
ax1.axvline(-0.2, color='gray', linestyle='--', linewidth=2)
ax1.axvline(0.2, color='gray', linestyle='--', linewidth=2)

ax2 = ax1.twinx()
sns.histplot(df_clusters_nhood, x="normalized_changed_cell_num", hue="type", kde=False, ax=ax2, stat="density")

plt.show()

from scipy.stats import norm

mu, std = norm.fit(df_clusters_nhood["normalized_changed_cell_num"])
weights = np.ones_like(df_clusters_nhood["normalized_changed_cell_num"]) / len(df_clusters_nhood["normalized_changed_cell_num"])

# df_clusters_nhood.pivot(columns="type", values="normalized_changed_cell_num").plot.hist(bins=100, weights=weights)
plt.hist(df_clusters_nhood["normalized_changed_cell_num"], bins=100, weights=weights)
# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)

plt.plot(x, p, 'k', linewidth=2)



plt.show()


df_clusters_nhood.pivot(columns="type", values="normalized_changed_cell_num").plot.hist(bins=100, normed=True)