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

GO_df_up_day0 = pd.read_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/GO_df_up_day0.csv")
GO_df_down_day0 = pd.read_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/GO_df_down_day0.csv")

GO_df_up_day0 = GO_df_up_day0[["term", "class", "cluster"]]
GO_df_up_day0 = GO_df_up_day0.drop_duplicates()

GO_df_down_day0 = GO_df_down_day0[["term", "class", "cluster"]]
GO_df_down_day0 = GO_df_down_day0.drop_duplicates()

## dep vs res

GO_df_up_day0_dep = GO_df_up_day0[GO_df_up_day0["cluster"].isin(dep_clusters_neuron_only)]
GO_df_up_day0_dep["cluster"] = GO_df_up_day0_dep["cluster"].astype(str)
GO_df_up_day0_dep_a = GO_df_up_day0_dep.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_day0_dep_b = GO_df_up_day0_dep.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_day0_dep = GO_df_up_day0_dep_a.merge(GO_df_up_day0_dep_b)
GO_df_up_day0_dep = GO_df_up_day0_dep.sort_values(by='cluster_counts', ascending=False)
GO_df_up_day0_dep.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_up_day0_dep.csv")

GO_df_down_day0_dep = GO_df_down_day0[GO_df_down_day0["cluster"].isin(dep_clusters_neuron_only)]
GO_df_down_day0_dep["cluster"] = GO_df_down_day0_dep["cluster"].astype(str)
GO_df_down_day0_dep_a = GO_df_down_day0_dep.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_day0_dep_b = GO_df_down_day0_dep.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_day0_dep = GO_df_down_day0_dep_a.merge(GO_df_down_day0_dep_b)
GO_df_down_day0_dep = GO_df_down_day0_dep.sort_values(by='cluster_counts', ascending=False)
GO_df_down_day0_dep.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_down_day0_dep.csv")

GO_df_up_day0_res = GO_df_up_day0[GO_df_up_day0["cluster"].isin(res_clusters_neuron_only)]
GO_df_up_day0_res["cluster"] = GO_df_up_day0_res["cluster"].astype(str)
GO_df_up_day0_res_a = GO_df_up_day0_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_day0_res_b = GO_df_up_day0_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_day0_res = GO_df_up_day0_res_a.merge(GO_df_up_day0_res_b)
GO_df_up_day0_res = GO_df_up_day0_res.sort_values(by='cluster_counts', ascending=False)
GO_df_up_day0_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_up_day0_res.csv")

GO_df_down_day0_res = GO_df_down_day0[GO_df_down_day0["cluster"].isin(res_clusters_neuron_only)]
GO_df_down_day0_res["cluster"] = GO_df_down_day0_res["cluster"].astype(str)
GO_df_down_day0_res_a = GO_df_down_day0_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_day0_res_b = GO_df_down_day0_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_day0_res = GO_df_down_day0_res_a.merge(GO_df_down_day0_res_b)
GO_df_down_day0_res = GO_df_down_day0_res.sort_values(by='cluster_counts', ascending=False)
GO_df_down_day0_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/dep_vs_res/GO_compare/GO_df_down_day0_res.csv")

comparison4 = {"depleted_up": set(GO_df_up_day0_dep.term.values),
               "depleted_down": set(GO_df_down_day0_dep.term.values),
               "resistant_up": set(GO_df_up_day0_res.term.values),
               "resistant_down": set(GO_df_down_day0_res.term.values)}

venn(comparison4)

## changed vs res

changed_clusters_neuron_only = dep_clusters_neuron_only + mov_clusters_neuron_only

GO_df_up_day0_changed = GO_df_up_day0[GO_df_up_day0["cluster"].isin(changed_clusters_neuron_only)]
GO_df_up_day0_changed["cluster"] = GO_df_up_day0_changed["cluster"].astype(str)
GO_df_up_day0_changed_a = GO_df_up_day0_changed.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_day0_changed_b = GO_df_up_day0_changed.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_day0_changed = GO_df_up_day0_changed_a.merge(GO_df_up_day0_changed_b)
GO_df_up_day0_changed = GO_df_up_day0_changed.sort_values(by='cluster_counts', ascending=False)
GO_df_up_day0_changed.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_up_day0_changed.csv")

GO_df_down_day0_changed = GO_df_down_day0[GO_df_down_day0["cluster"].isin(changed_clusters_neuron_only)]
GO_df_down_day0_changed["cluster"] = GO_df_down_day0_changed["cluster"].astype(str)
GO_df_down_day0_changed_a = GO_df_down_day0_changed.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_day0_changed_b = GO_df_down_day0_changed.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_day0_changed = GO_df_down_day0_changed_a.merge(GO_df_down_day0_changed_b)
GO_df_down_day0_changed = GO_df_down_day0_changed.sort_values(by='cluster_counts', ascending=False)
GO_df_down_day0_changed.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_down_day0_changed.csv")

GO_df_up_day0_res = GO_df_up_day0[GO_df_up_day0["cluster"].isin(res_clusters_neuron_only)]
GO_df_up_day0_res["cluster"] = GO_df_up_day0_res["cluster"].astype(str)
GO_df_up_day0_res_a = GO_df_up_day0_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_day0_res_b = GO_df_up_day0_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_day0_res = GO_df_up_day0_res_a.merge(GO_df_up_day0_res_b)
GO_df_up_day0_res = GO_df_up_day0_res.sort_values(by='cluster_counts', ascending=False)
GO_df_up_day0_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_up_day0_res.csv")

GO_df_down_day0_res = GO_df_down_day0[GO_df_down_day0["cluster"].isin(res_clusters_neuron_only)]
GO_df_down_day0_res["cluster"] = GO_df_down_day0_res["cluster"].astype(str)
GO_df_down_day0_res_a = GO_df_down_day0_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_day0_res_b = GO_df_down_day0_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_day0_res = GO_df_down_day0_res_a.merge(GO_df_down_day0_res_b)
GO_df_down_day0_res = GO_df_down_day0_res.sort_values(by='cluster_counts', ascending=False)
GO_df_down_day0_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/changed_vs_res/GO_compare/GO_df_down_day0_res.csv")

comparison4 = {"changed_up": set(GO_df_up_day0_changed.term.values),
               "changed_down": set(GO_df_down_day0_changed.term.values),
               "resistant_up": set(GO_df_up_day0_res.term.values),
               "resistant_down": set(GO_df_down_day0_res.term.values)}

venn(comparison4)

## mov vs res

GO_df_up_day0_mov = GO_df_up_day0[GO_df_up_day0["cluster"].isin(mov_clusters_neuron_only)]
GO_df_up_day0_mov["cluster"] = GO_df_up_day0_mov["cluster"].astype(str)
GO_df_up_day0_mov_a = GO_df_up_day0_mov.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_day0_mov_b = GO_df_up_day0_mov.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_day0_mov = GO_df_up_day0_mov_a.merge(GO_df_up_day0_mov_b)
GO_df_up_day0_mov = GO_df_up_day0_mov.sort_values(by='cluster_counts', ascending=False)
GO_df_up_day0_mov.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_up_day0_mov.csv")

GO_df_down_day0_mov = GO_df_down_day0[GO_df_down_day0["cluster"].isin(mov_clusters_neuron_only)]
GO_df_down_day0_mov["cluster"] = GO_df_down_day0_mov["cluster"].astype(str)
GO_df_down_day0_mov_a = GO_df_down_day0_mov.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_day0_mov_b = GO_df_down_day0_mov.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_day0_mov = GO_df_down_day0_mov_a.merge(GO_df_down_day0_mov_b)
GO_df_down_day0_mov = GO_df_down_day0_mov.sort_values(by='cluster_counts', ascending=False)
GO_df_down_day0_mov.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_down_day0_mov.csv")

GO_df_up_day0_res = GO_df_up_day0[GO_df_up_day0["cluster"].isin(res_clusters_neuron_only)]
GO_df_up_day0_res["cluster"] = GO_df_up_day0_res["cluster"].astype(str)
GO_df_up_day0_res_a = GO_df_up_day0_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_up_day0_res_b = GO_df_up_day0_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_up_day0_res = GO_df_up_day0_res_a.merge(GO_df_up_day0_res_b)
GO_df_up_day0_res = GO_df_up_day0_res.sort_values(by='cluster_counts', ascending=False)
GO_df_up_day0_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_up_day0_res.csv")

GO_df_down_day0_res = GO_df_down_day0[GO_df_down_day0["cluster"].isin(res_clusters_neuron_only)]
GO_df_down_day0_res["cluster"] = GO_df_down_day0_res["cluster"].astype(str)
GO_df_down_day0_res_a = GO_df_down_day0_res.groupby(['term', 'class'])['cluster'].agg(lambda x: ', '.join(x)).reset_index(name='cluster_list')
GO_df_down_day0_res_b = GO_df_down_day0_res.groupby(['term', 'class']).size().reset_index(name='cluster_counts')
GO_df_down_day0_res = GO_df_down_day0_res_a.merge(GO_df_down_day0_res_b)
GO_df_down_day0_res = GO_df_down_day0_res.sort_values(by='cluster_counts', ascending=False)
GO_df_down_day0_res.to_csv(r"/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8_corrected/mov_vs_res/GO_compare/GO_df_down_day0_res.csv")

comparison4 = {"moved_up": set(GO_df_up_day0_mov.term.values),
               "moved_down": set(GO_df_down_day0_mov.term.values),
               "resistant_up": set(GO_df_up_day0_res.term.values),
               "resistant_down": set(GO_df_down_day0_res.term.values)}

venn(comparison4)
