# bare-bones
Bare-bones Gene Ontology enrichment is written in Python3 to conduct hyper- geometric testing between experimental *Arabidopsis thaliana* gene sets and *A. thaliana*, GOSLIM Gene Ontology database. Bare-bones has been tested on a Windows 10 machine. 

Use `example.py` as a guide on how to implement `structures.py`, `set_parameters.py` and `enrichment.py`.

| Class | Description |
| --- | --- |
| `structures.py`     | Initializes data structures for significant enrichment test(s). |
| `set_parameters.py` | Helper class for passing parameters from Gene Ontology data structure. |
| `enrichment.py`     | Conducts enrichment test and exports enrichment results. |

Gene set enrichment results are written to a file, following the structure below. Additionally, any unannotated genes in input gene sets are output to files corresponding to the gene set. 

Below is an example of an unsorted enrichment result; NA is listed in the indentifiers column if there is no significant enrichment as determined by the significance cut- off provided.

|     aspect    |     go_term   |  significance |  identifiers  |
| ------------- | ------------- | ------------- | ------------- |
| P             | GO:0006355    | 4.343788313513017e-54 | AT1G01010 AT1G01030 AT1G01060 ... |
| P             | GO:0032366    | 0.13531096220565436 | NA  |
| C             | GO:0000123    | 0.02967922876604696| AT1G02740 |
| C             | GO:0001401    | 0.4095373894550747 | NA  |
| F             | GO:0003725    | 0.014840807240371415 | AT1G01040 |
| F             | GO:0004525    | 0.051216165118646684 | NA |
