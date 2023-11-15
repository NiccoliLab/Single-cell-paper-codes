import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from venn import venn
from matplotlib_venn import venn3

dep_clusters = pd.read_csv("/home/sdxucl/Desktop/dep_clusters.csv")
dep_clusters = dep_clusters["clusters"].to_list()

res_clusters = pd.read_csv("/home/sdxucl/Desktop/res_clusters.csv")
res_clusters = res_clusters["clusters"].to_list()

mov_clusters = pd.read_csv("/home/sdxucl/Desktop/mov_clusters.csv")
mov_clusters = mov_clusters["clusters"].to_list()

glia_clusters = set({72, 48, 74, 91, 137, 107})

dep_clusters_neuron_only = list(set(dep_clusters) - glia_clusters)
res_clusters_neuron_only = list(set(res_clusters) - glia_clusters)
mov_clusters_neuron_only = list(set(mov_clusters) - glia_clusters)

GO_df_up_DE = pd.read_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/GO_df_up_DE.csv")
GO_df_down_DE = pd.read_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/GO_df_down_DE.csv")

GO_df_up_DE = GO_df_up_DE[["term", "class", "cluster"]]
GO_df_up_DE = GO_df_up_DE.drop_duplicates()

GO_df_down_DE = GO_df_down_DE[["term", "class", "cluster"]]
GO_df_down_DE = GO_df_down_DE.drop_duplicates()

## dep vs res

GO_df_up_DE_dep = GO_df_up_DE[GO_df_up_DE["cluster"].isin(dep_clusters_neuron_only)]
GO_df_up_DE_dep["cluster"] = GO_df_up_DE_dep["cluster"].astype(str)
GO_df_up_DE_dep_a = GO_df_up_DE_dep.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_DE_dep_b = GO_df_up_DE_dep.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_DE_dep = GO_df_up_DE_dep_a.merge(GO_df_up_DE_dep_b)
GO_df_up_DE_dep = GO_df_up_DE_dep.sort_values(by='cluster_counts', ascending=False)
GO_df_up_DE_dep.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_up_DE_dep.csv")

GO_df_down_DE_dep = GO_df_down_DE[GO_df_down_DE["cluster"].isin(dep_clusters_neuron_only)]
GO_df_down_DE_dep["cluster"] = GO_df_down_DE_dep["cluster"].astype(str)
GO_df_down_DE_dep_a = GO_df_down_DE_dep.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_DE_dep_b = GO_df_down_DE_dep.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_DE_dep = GO_df_down_DE_dep_a.merge(GO_df_down_DE_dep_b)
GO_df_down_DE_dep = GO_df_down_DE_dep.sort_values(by='cluster_counts', ascending=False)
GO_df_down_DE_dep.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_down_DE_dep.csv")

GO_df_up_DE_res = GO_df_up_DE[GO_df_up_DE["cluster"].isin(res_clusters_neuron_only)]
GO_df_up_DE_res["cluster"] = GO_df_up_DE_res["cluster"].astype(str)
GO_df_up_DE_res_a = GO_df_up_DE_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_DE_res_b = GO_df_up_DE_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_DE_res = GO_df_up_DE_res_a.merge(GO_df_up_DE_res_b)
GO_df_up_DE_res = GO_df_up_DE_res.sort_values(by='cluster_counts', ascending=False)
GO_df_up_DE_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_up_DE_res.csv")

GO_df_down_DE_res = GO_df_down_DE[GO_df_down_DE["cluster"].isin(res_clusters_neuron_only)]
GO_df_down_DE_res["cluster"] = GO_df_down_DE_res["cluster"].astype(str)
GO_df_down_DE_res_a = GO_df_down_DE_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_DE_res_b = GO_df_down_DE_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_DE_res = GO_df_down_DE_res_a.merge(GO_df_down_DE_res_b)
GO_df_down_DE_res = GO_df_down_DE_res.sort_values(by='cluster_counts', ascending=False)
GO_df_down_DE_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_down_DE_res.csv")

comparison4 = {"depleted_up": set(GO_df_up_DE_dep.term.values),
               "depleted_down": set(GO_df_down_DE_dep.term.values),
               "resistant_up": set(GO_df_up_DE_res.term.values),
               "resistant_down": set(GO_df_down_DE_res.term.values)}

venn(comparison4)

## changed vs res

changed_clusters_neuron_only = dep_clusters_neuron_only + mov_clusters_neuron_only

GO_df_up_DE_changed = GO_df_up_DE[GO_df_up_DE["cluster"].isin(changed_clusters_neuron_only)]
GO_df_up_DE_changed["cluster"] = GO_df_up_DE_changed["cluster"].astype(str)
GO_df_up_DE_changed_a = GO_df_up_DE_changed.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_DE_changed_b = GO_df_up_DE_changed.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_DE_changed = GO_df_up_DE_changed_a.merge(GO_df_up_DE_changed_b)
GO_df_up_DE_changed = GO_df_up_DE_changed.sort_values(by='cluster_counts', ascending=False)
GO_df_up_DE_changed.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_up_DE_changed.csv")

GO_df_down_DE_changed = GO_df_down_DE[GO_df_down_DE["cluster"].isin(changed_clusters_neuron_only)]
GO_df_down_DE_changed["cluster"] = GO_df_down_DE_changed["cluster"].astype(str)
GO_df_down_DE_changed_a = GO_df_down_DE_changed.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_DE_changed_b = GO_df_down_DE_changed.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_DE_changed = GO_df_down_DE_changed_a.merge(GO_df_down_DE_changed_b)
GO_df_down_DE_changed = GO_df_down_DE_changed.sort_values(by='cluster_counts', ascending=False)
GO_df_down_DE_changed.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_down_DE_changed.csv")

GO_df_up_DE_res = GO_df_up_DE[GO_df_up_DE["cluster"].isin(res_clusters_neuron_only)]
GO_df_up_DE_res["cluster"] = GO_df_up_DE_res["cluster"].astype(str)
GO_df_up_DE_res_a = GO_df_up_DE_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_DE_res_b = GO_df_up_DE_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_DE_res = GO_df_up_DE_res_a.merge(GO_df_up_DE_res_b)
GO_df_up_DE_res = GO_df_up_DE_res.sort_values(by='cluster_counts', ascending=False)
GO_df_up_DE_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_up_DE_res.csv")

GO_df_down_DE_res = GO_df_down_DE[GO_df_down_DE["cluster"].isin(res_clusters_neuron_only)]
GO_df_down_DE_res["cluster"] = GO_df_down_DE_res["cluster"].astype(str)
GO_df_down_DE_res_a = GO_df_down_DE_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_DE_res_b = GO_df_down_DE_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_DE_res = GO_df_down_DE_res_a.merge(GO_df_down_DE_res_b)
GO_df_down_DE_res = GO_df_down_DE_res.sort_values(by='cluster_counts', ascending=False)
GO_df_down_DE_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_down_DE_res.csv")

comparison4 = {"changed_up": set(GO_df_up_DE_changed.term.values),
               "changed_down": set(GO_df_down_DE_changed.term.values),
               "resistant_up": set(GO_df_up_DE_res.term.values),
               "resistant_down": set(GO_df_down_DE_res.term.values)}

venn(comparison4)

## mov vs res

GO_df_up_DE_mov = GO_df_up_DE[GO_df_up_DE["cluster"].isin(mov_clusters_neuron_only)]
GO_df_up_DE_mov["cluster"] = GO_df_up_DE_mov["cluster"].astype(str)
GO_df_up_DE_mov_a = GO_df_up_DE_mov.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_DE_mov_b = GO_df_up_DE_mov.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_DE_mov = GO_df_up_DE_mov_a.merge(GO_df_up_DE_mov_b)
GO_df_up_DE_mov = GO_df_up_DE_mov.sort_values(by='cluster_counts', ascending=False)
GO_df_up_DE_mov.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_up_DE_mov.csv")

GO_df_down_DE_mov = GO_df_down_DE[GO_df_down_DE["cluster"].isin(mov_clusters_neuron_only)]
GO_df_down_DE_mov["cluster"] = GO_df_down_DE_mov["cluster"].astype(str)
GO_df_down_DE_mov_a = GO_df_down_DE_mov.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_DE_mov_b = GO_df_down_DE_mov.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_DE_mov = GO_df_down_DE_mov_a.merge(GO_df_down_DE_mov_b)
GO_df_down_DE_mov = GO_df_down_DE_mov.sort_values(by='cluster_counts', ascending=False)
GO_df_down_DE_mov.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_down_DE_mov.csv")

GO_df_up_DE_res = GO_df_up_DE[GO_df_up_DE["cluster"].isin(res_clusters_neuron_only)]
GO_df_up_DE_res["cluster"] = GO_df_up_DE_res["cluster"].astype(str)
GO_df_up_DE_res_a = GO_df_up_DE_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_DE_res_b = GO_df_up_DE_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_DE_res = GO_df_up_DE_res_a.merge(GO_df_up_DE_res_b)
GO_df_up_DE_res = GO_df_up_DE_res.sort_values(by='cluster_counts', ascending=False)
GO_df_up_DE_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_up_DE_res.csv")

GO_df_down_DE_res = GO_df_down_DE[GO_df_down_DE["cluster"].isin(res_clusters_neuron_only)]
GO_df_down_DE_res["cluster"] = GO_df_down_DE_res["cluster"].astype(str)
GO_df_down_DE_res_a = GO_df_down_DE_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_DE_res_b = GO_df_down_DE_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_DE_res = GO_df_down_DE_res_a.merge(GO_df_down_DE_res_b)
GO_df_down_DE_res = GO_df_down_DE_res.sort_values(by='cluster_counts', ascending=False)
GO_df_down_DE_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_down_DE_res.csv")

comparison4 = {"moved_up": set(GO_df_up_DE_mov.term.values),
               "moved_down": set(GO_df_down_DE_mov.term.values),
               "resistant_up": set(GO_df_up_DE_res.term.values),
               "resistant_down": set(GO_df_down_DE_res.term.values)}

venn(comparison4)

plt.rcParams.update({'font.size': 22})
ax = down.plot.bar(x="term", y="cluster_counts", color=[colors[i] for i in down["class"]])
plt.xticks(rotation=45, ha='right')
plt.ylabel("cluster_counts")
red = mpatches.Patch(color='red', label='biological_process')
green = mpatches.Patch(color='green', label='cellular_component')
blue = mpatches.Patch(color='blue', label="molecular_function")
plt.legend(handles=[red, green, blue])
