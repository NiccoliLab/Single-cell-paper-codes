import numpy as np
import pandas as pd

fly_human = pd.read_csv("/home/sdxucl/Desktop/fly_human_ortholog.txt")

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

human_genes = pd.read_csv("/home/sdxucl/Desktop/human_oMN_Onuf.csv")

human_OMN_up_only = human_genes[["OMN only \nupregulated"]].dropna()
human_OMN_up_only = human_OMN_up_only.rename(columns={"OMN only \nupregulated": "Human gene name"})
human_OMN_down_only = human_genes[["OMN only \ndownregulated"]].dropna()
human_OMN_down_only = human_OMN_down_only.rename(columns={"OMN only \ndownregulated": "Human gene name"})

human_Onuf_up_only = human_genes[["Onuf only \nupregulated"]].dropna()
human_Onuf_up_only = human_Onuf_up_only.rename(columns={"Onuf only \nupregulated": "Human gene name"})
human_Onuf_down_only = human_genes[["Onuf only \ndownregulated"]].dropna()
human_Onuf_down_only = human_Onuf_down_only.rename(columns={"Onuf only \ndownregulated": "Human gene name"})

human_commonly_up = human_genes[["Commonly \nupregulated"]].dropna()
human_commonly_up = human_commonly_up.rename(columns={"Commonly \nupregulated": "Human gene name"})
human_commonly_down = human_genes[["Commonly \ndownregulated"]].dropna()
human_commonly_down = human_commonly_down.rename(columns={"Commonly \ndownregulated": "Human gene name"})

human_OMN_up = pd.concat([human_OMN_up_only, human_commonly_up])
human_OMN_down = pd.concat([human_OMN_down_only, human_commonly_down])

human_Onuf_up = pd.concat([human_Onuf_up_only, human_commonly_up])
human_Onuf_down = pd.concat([human_Onuf_down_only, human_commonly_down])

human_OMN_up_FLY = human_OMN_up_only.merge(fly_human, how="left", on="Human gene name").dropna()
human_OMN_down_FLY = human_OMN_down_only.merge(fly_human, how="left", on="Human gene name").dropna()
human_Onuf_up_FLY = human_Onuf_up_only.merge(fly_human, how="left", on="Human gene name").dropna()
human_Onuf_down_FLY = human_Onuf_down.merge(fly_human, how="left", on="Human gene name").dropna()
human_commonly_up_FLY = human_commonly_up.merge(fly_human, how="left", on="Human gene name").dropna()
human_commonly_down_FLY = human_commonly_down.merge(fly_human, how="left", on="Human gene name").dropna()

day0_up_res = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/up/DE_genes_up_2.csv")
for x in res_clusters_neuron_only:
    day0_up_res_x = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/up/DE_genes_up_" + str(x) + ".csv")
    day0_up_res = pd.concat([day0_up_res, day0_up_res_x])
day0_up_res_morethan5 = day0_up_res.groupby("names").count()
day0_up_res_morethan5 = day0_up_res_morethan5[["group"]]
day0_up_res_morethan5 = day0_up_res_morethan5.loc[day0_up_res_morethan5["group"] >= 5]

day0_up_dep = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/up/DE_genes_up_0.csv")
for x in dep_clusters_neuron_only:
    day0_up_dep_x = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/up/DE_genes_up_" + str(x) + ".csv")
    day0_up_dep = pd.concat([day0_up_dep, day0_up_dep_x])
day0_up_dep_morethan5 = day0_up_dep.groupby("names").count()
day0_up_dep_morethan5 = day0_up_dep_morethan5[["group"]]
day0_up_dep_morethan5 = day0_up_dep_morethan5.loc[day0_up_dep_morethan5["group"] >= 5]

day0_down_res = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/down/DE_genes_down_2.csv")
for x in res_clusters_neuron_only:
    day0_down_res_x = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/down/DE_genes_down_" + str(x) + ".csv")
    day0_down_res = pd.concat([day0_down_res, day0_down_res_x])

day0_down_dep = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/down/DE_genes_down_0.csv")
for x in dep_clusters_neuron_only:
    day0_down_dep_x = pd.read_csv("/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/down/DE_genes_down_" + str(x) + ".csv")
    day0_down_dep = pd.concat([day0_down_dep, day0_down_dep_x])

# day0_up_res_unique = set(day0_up_res["names"]) - (set(day0_up_res["names"]) & set(day0_up_dep["names"]))
# day0_up_dep_unique = set(day0_up_dep["names"]) - (set(day0_up_res["names"]) & set(day0_up_dep["names"]))

day0_up_res_unique = set(day0_up_res_morethan5.index) - (set(day0_up_res_morethan5.index) & set(day0_up_dep_morethan5.index))
day0_up_dep_unique = set(day0_up_dep_morethan5.index) - (set(day0_up_res_morethan5.index) & set(day0_up_dep_morethan5.index))

day0_down_res_unique = set(day0_down_res["names"]) - (set(day0_down_res["names"]) & set(day0_down_dep["names"]))
day0_down_dep_unique = set(day0_down_dep["names"]) - (set(day0_down_res["names"]) & set(day0_down_dep["names"]))

from matplotlib_venn import venn2
from matplotlib_venn import venn3

out = venn2([set(human_OMN_up["Human gene name"]), set(human_Onuf_up["Human gene name"])], ('Human_OMN_up', 'Human_Onuf_up'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn2([set(human_OMN_down["Human gene name"]), set(human_Onuf_down["Human gene name"])], ('Human_OMN_down', 'Human_Onuf_down'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_OMN_up_FLY["Gene name"]), set(day0_up_res["names"]), set(day0_up_dep["names"])], ('Human_OMN_up_FLY', 'fly_res_day0_high', 'fly_dep_day0_high'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_Onuf_up_FLY["Gene name"]), set(day0_up_res["names"]), set(day0_up_dep["names"])], ('Human_Onuf_up_FLY', 'fly_res_day0_high', 'fly_dep_day0_high'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_commonly_up_FLY["Gene name"]), set(day0_up_res["names"]), set(day0_up_dep["names"])], ('Human_commonly_up_FLY', 'fly_res_day0_high', 'fly_dep_day0_high'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_OMN_up_FLY["Gene name"]), day0_up_res_unique, day0_up_dep_unique], ('Human_OMN_up_FLY', 'fly_res_day0_high_unique_over5', 'fly_dep_day0_high_unique_over5'))
for text in out.set_labels:
    text.set_fontsize(10)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(10)

overlap_OMN_up_day0_resistant = pd.DataFrame(set(human_OMN_up_FLY["Gene name"]) & day0_up_res_unique)
overlap_OMN_up_day0_resistant = overlap_OMN_up_day0_resistant.rename(columns={0: "Gene name"})
overlap_OMN_up_day0_resistant_with_Human = overlap_OMN_up_day0_resistant.merge(human_OMN_up_FLY, how="left", on="Gene name")
overlap_OMN_up_day0_resistant.to_csv("/home/sdxucl/Desktop/Teresa_Chicago_conference/Comparison_oMNs_Onufs_human/overlap_OMN_up_day0_resistant_over5.csv")
overlap_OMN_up_day0_resistant_with_Human.to_csv("/home/sdxucl/Desktop/Teresa_Chicago_conference/Comparison_oMNs_Onufs_human/overlap_OMN_up_day0_resistant_over5_with_Human.csv")

out = venn3([set(human_Onuf_up_FLY["Gene name"]), day0_up_res_unique, day0_up_dep_unique], ('Human_Onuf_up_FLY', 'fly_res_day0_high_unique', 'fly_dep_day0_high_unique'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

overlap_Onuf_up_day0_resistant = pd.DataFrame(set(human_Onuf_up_FLY["Gene name"]) & day0_up_res_unique)
overlap_Onuf_up_day0_resistant = overlap_Onuf_up_day0_resistant.rename(columns={0: "Gene name"})
overlap_Onuf_up_day0_resistant_with_Human = overlap_Onuf_up_day0_resistant.merge(human_Onuf_up_FLY, how="left", on="Gene name")
overlap_Onuf_up_day0_resistant.to_csv("/home/sdxucl/Desktop/Teresa_Chicago_conference/Comparison_oMNs_Onufs_human/overlap_Onuf_up_day0_resistant_over5.csv")
overlap_Onuf_up_day0_resistant_with_Human.to_csv("/home/sdxucl/Desktop/Teresa_Chicago_conference/Comparison_oMNs_Onufs_human/overlap_Onuf_up_day0_resistant_over5_with_Human.csv")

out = venn3([set(human_commonly_up_FLY["Gene name"]), day0_up_res_unique, day0_up_dep_unique], ('Human_commonly_up_FLY', 'fly_res_day0_high_unique', 'fly_dep_day0_high_unique'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

overlap_commonly_up_day0_resistant = pd.DataFrame(set(human_commonly_up_FLY["Gene name"]) & day0_up_res_unique)
overlap_commonly_up_day0_resistant.to_csv("/home/sdxucl/Desktop/Teresa_Chicago_conference/Comparison_oMNs_Onufs_human/overlap_commonly_up_day0_resistant_over5.csv")

out = venn3([set(human_OMN_down_FLY["Gene name"]), set(day0_down_res["names"]), set(day0_down_dep["names"])], ('Human_OMN_down_FLY', 'fly_res_day0_high', 'fly_dep_day0_high'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_Onuf_down_FLY["Gene name"]), set(day0_down_res["names"]), set(day0_down_dep["names"])], ('Human_Onuf_down_FLY', 'fly_res_day0_high', 'fly_dep_day0_high'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_commonly_down_FLY["Gene name"]), set(day0_down_res["names"]), set(day0_down_dep["names"])], ('Human_commonly_down_FLY', 'fly_res_day0_high', 'fly_dep_day0_high'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_OMN_down_FLY["Gene name"]), day0_down_res_unique, day0_down_dep_unique], ('Human_OMN_down_FLY', 'fly_res_day0_high_unique', 'fly_dep_day0_high_unique'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_Onuf_down_FLY["Gene name"]), day0_down_res_unique, day0_down_dep_unique], ('Human_Onuf_down_FLY', 'fly_res_day0_low_unique', 'fly_dep_day0_low_unique'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

out = venn3([set(human_commonly_down_FLY["Gene name"]), day0_down_res_unique, day0_down_dep_unique], ('Human_commonly_down_FLY', 'fly_res_day0_low_unique', 'fly_dep_day0_low_unique'))
for text in out.set_labels:
    text.set_fontsize(20)
for x in range(len(out.subset_labels)):
    if out.subset_labels[x] is not None:
        out.subset_labels[x].set_fontsize(20)

overlap_commonly_down_day0_resistant = pd.DataFrame(set(human_commonly_up_FLY["Gene name"]) & day0_down_res_unique)
overlap_commonly_down_day0_resistant.to_csv("/home/sdxucl/Desktop/Teresa_Chicago_conference/Comparison_oMNs_Onufs_human/overlap_commonly_down_day0_resistant.csv")

