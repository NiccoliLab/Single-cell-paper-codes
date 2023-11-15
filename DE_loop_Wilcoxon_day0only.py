import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=150, facecolor='white')

corr_data = sc.read('c9_aggr_mnn_correction.h5ad')
adata_filter = sc.read('c9_aggr_prelim_filter.h5ad')
sc.pp.normalize_total(adata_filter, target_sum=1e4)
sc.pp.log1p(adata_filter)
corr_data.raw = adata_filter

sc.pp.neighbors(corr_data, n_neighbors=20, n_pcs=40)
sc.tl.leiden(corr_data, key_added='clusters', resolution=10.0)

sc.pl.umap(corr_data, color='clusters', s=50, frameon=False, palette=color_list, legend_fontsize="small")
# refer to "Generating_random_colors" for color_list

## Create a dataframe counting cell percentage along the time points in each cluster
corr_data.obs.groupby("day").apply(len)

norm_cellnum = pd.DataFrame(columns=['0day', '2day', '4day', '6day', '8day'])
norm_cellnum['0day'] = corr_data.obs[corr_data.obs['day'].isin(['0'])].groupby("clusters").apply(len) / \
                       (corr_data.obs["day"] == "0").sum()
norm_cellnum['2day'] = corr_data.obs[corr_data.obs['day'].isin(['2'])].groupby("clusters").apply(len) / \
                       (corr_data.obs["day"] == "2").sum()
norm_cellnum['4day'] = corr_data.obs[corr_data.obs['day'].isin(['4'])].groupby("clusters").apply(len) / \
                       (corr_data.obs["day"] == "4").sum()
norm_cellnum['6day'] = corr_data.obs[corr_data.obs['day'].isin(['6'])].groupby("clusters").apply(len) / \
                       (corr_data.obs["day"] == "6").sum()
norm_cellnum['8day'] = corr_data.obs[corr_data.obs['day'].isin(['8'])].groupby("clusters").apply(len) / \
                       (corr_data.obs["day"] == "8").sum()

corr_data_day0 = corr_data[corr_data.obs["day"] == "0"]  # extract cells at day 0

day0_clusters = corr_data_day0.obs.groupby("clusters").apply(len)

corr_data_day0 = corr_data_day0[corr_data_day0.obs["clusters"] != "38"]  # exclude clusters containing only 1 cell
corr_data_day0 = corr_data_day0[corr_data_day0.obs["clusters"] != "46"]  # exclude clusters containing only 1 cell

sc.tl.rank_genes_groups(corr_data_day0, groupby="clusters", method="wilcoxon", use_raw=True)

DE_genes = sc.get.rank_genes_groups_df(corr_data_day0, group=None)
DE_genes = DE_genes.drop(DE_genes[DE_genes.pvals_adj >= 0.05].index)
DE_genes.to_csv(r"/home/sdxucl/Desktop/DE_day0_clusters_WC.csv")

for x in list(range(0, 38)) + list(range(39, 46)) + list(range(47, 213)):
    DE_genes_up = DE_genes[(DE_genes["group"] == str(x)) & (DE_genes["logfoldchanges"] > 0)]
    DE_genes_down = DE_genes[(DE_genes["group"] == str(x)) & (DE_genes["logfoldchanges"] < 0)]
    DE_genes_up.to_csv("/home/sdxucl/Desktop/DE_day0_WC/up/DE_genes_up_" + str(x) + ".csv")
    DE_genes_down.to_csv("/home/sdxucl/Desktop/DE_day0_WC/down/DE_genes_down_" + str(x) + ".csv")


import numpy as np
import pandas as pd
import scanpy as sc

from goatools.base import download_go_basic_obo

obo_fname = download_go_basic_obo()

from goatools.base import download_ncbi_associations

fin_gene2go = download_ncbi_associations()

from goatools.obo_parser import GODag

obodag = GODag("go-basic.obo")

from __future__ import print_function
from goatools.anno.genetogo_reader import Gene2GoReader

objanno = Gene2GoReader(fin_gene2go, taxids=[7227])
ns2assoc = objanno.get_ns2assc()

for nspc, id2gos in ns2assoc.items():
    print("{NS} {N:,} annotated fly genes".format(NS=nspc, N=len(id2gos)))

from genes_ncbi_7227_proteincoding import GENEID2NT as GeneID2nt_dro

print(len(GeneID2nt_dro))

from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

goeaobj = GOEnrichmentStudyNS(
    GeneID2nt_dro.keys(),
    ns2assoc,
    obodag,
    propagate_counts=False,
    alpha=0.05,
    methods=['fdr_bh']
)

mapper = {}

for key in GeneID2nt_dro:
    mapper[GeneID2nt_dro[key].Symbol] = GeneID2nt_dro[key].GeneID

inv_map = {v: k for k, v in mapper.items()}

## mapper = pd.DataFrame.from_dict(mapper, orient='index')
## TT_44_8vs0_down_geneID = mapper[mapper[0].isin(TT_44_8vs0_down["names"].tolist())]
## TT_44_8vs0_down_input = dict(zip(TT_44_8vs0_down_geneID.index.tolist(), TT_44_8vs0_down_geneID[0].tolist())).keys()

## goea_results_all = goeaobj.run_study(TT_44_8vs0_down_input)
## goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

GO_items = []

temp = goeaobj.ns2objgoea['BP'].assoc
for item in temp:
    GO_items += temp[item]

temp = goeaobj.ns2objgoea['CC'].assoc
for item in temp:
    GO_items += temp[item]

temp = goeaobj.ns2objgoea['MF'].assoc
for item in temp:
    GO_items += temp[item]


def go_it(test_genes):
    print(f'input genes: {len(test_genes)}')

    mapped_genes = []
    for gene in test_genes:
        try:
            mapped_genes.append(mapper[gene])
        except:
            pass
    print(f'mapped genes: {len(mapped_genes)}')

    goea_results_all = goeaobj.run_study(mapped_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh, \
                                          x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO),
                                          list(map(lambda y: inv_map[y], x.study_items)), \
                                          ], goea_results_sig)),
                      columns=['GO', 'term', 'class', 'p', 'p_corr', 'n_genes', \
                               'n_study', 'n_go', 'study_genes'])

    GO = GO[GO.n_genes > 1]
    return GO

for x in list(range(0, 38)) + list(range(39, 46)) + list(range(47, 213)):
    input_up = pd.read_csv(
        "/home/sdxucl/Desktop/DE_day0_WC/up/DE_genes_up_" + str(x) + ".csv")
    input_down = pd.read_csv(
        "/home/sdxucl/Desktop/DE_day0_WC/down/DE_genes_down_" + str(x) + ".csv")
    df_up = go_it(input_up.names.values)
    df_up['per'] = df_up.n_genes / df_up.n_go
    df_up.to_csv("/home/sdxucl/Desktop/DE_day0_WC/GO/up/GO_up_" + str(x) + ".csv")
    df_down = go_it(input_down.names.values)
    df_down['per'] = df_down.n_genes / df_down.n_go
    df_down.to_csv("/home/sdxucl/Desktop/DE_day0_WC/GO/down/GO_down_" + str(x) + ".csv")
