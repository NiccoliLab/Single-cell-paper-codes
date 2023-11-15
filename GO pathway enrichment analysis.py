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


day = [2, 4, 6, 8]

for x in range(189, 204): # 82, 188, 204, 207
    for y in day:
        input_up = pd.read_csv(
            "/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/DE/DE_genes_up_" + str(x) + "_" + str(y) + "vs0.csv")
        input_down = pd.read_csv(
            "/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/DE/DE_genes_down_" + str(x) + "_" + str(y) + "vs0.csv")
        df_up = go_it(input_up.names.values)
        df_up['per'] = df_up.n_genes / df_up.n_go
        df_up.to_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8/DE/GO_up_" + str(x) + "_" + str(y) + "vs0.csv")
        df_down = go_it(input_down.names.values)
        df_down['per'] = df_down.n_genes / df_down.n_go
        df_down.to_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8/DE/GO_down_" + str(x) + "_" + str(y) + "vs0.csv")

for x in range(0, 213):
    for y in day:
        input_up = pd.read_csv(
            "/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/up/DE_genes_up_" + str(x) + ".csv")
        input_down = pd.read_csv(
            "/home/sdxucl/Desktop/DE_loop_wilcoxon_ex_b1d8/day0/down/DE_genes_down_" + str(x) + ".csv")
        df_up = go_it(input_up.names.values)
        df_up['per'] = df_up.n_genes / df_up.n_go
        df_up.to_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8/day0/GO_up_" + str(x) + ".csv")
        df_down = go_it(input_down.names.values)
        df_down['per'] = df_down.n_genes / df_down.n_go
        df_down.to_csv("/home/sdxucl/Desktop/GO_string_analysis_ex_b1d8/day0/GO_down_" + str(x) + ".csv")
