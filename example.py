from structures import Struct
from enrichment import G
# link to Gene Ontology annotations.
link = 'https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt'
# place file name(s) containing genes for testing.
test_list = ['gene_list.tsv']
# implement the Struct class.
struct = Struct(t = test_list, l = link)
# conduct Gene Ontology enrichment with significance cut- off.
go = G(struct, 0.05)
