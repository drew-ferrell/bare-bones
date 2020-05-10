from structures import Struct
from enrichment import G
link = 'https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt'
test_list = ['gene_list.tsv']
struct = Struct(t = test_list, l = link)
go = G(struct, 0.05)
