from set_parameters import Param
from errno import ENOENT
from requests import get
from os import getcwd, strerror
from os.path import basename, exists
from numpy import loadtxt, str
from csv import writer

class Struct():
    def __init__(self, t=None, l=None):
        '''
        t: file paths list
        l: Gene Ontology database link
        '''
        self.test_list = t
        self.link = l
        self.go_table = dict()

        self.database_genes = set()
        self.unannotated_genes = list()

        self.set_go_table()
        self.get_test_list()
        self.export_unannotated_genes()
               
    def get_test_list(self):
        '''
        load gene list contents to memory
        '''
        non_empty_list = self.non_empty_list()
        if non_empty_list:
            self.test_list = self.set_test_cont()
        else:
            raise Warning('Unexpected data type, or empty list:', self.test_list)
        return 
    
    def non_empty_list(self):
        '''
        generalizable data- type testing
        '''
        return self.non_empty(list)
        
    def non_empty(self, d):
        '''
        determines test list structure

        d: data type
        '''
        test_list = self.test_list
        empty = len(test_list) == 0
        is_data_type = (type(test_list) == d)
        if is_data_type and empty:
            raise Warning('Unexpected number of gene lists.')
        if is_data_type and not empty:
            return True
        return 
    
    def set_test_cont(self):
        '''
        insert genes into hash table as values
        '''
        test_list = self.test_list
        # construct gene set hash table
        construct_hash_table = self.set_test_hash
        test_structure = {test_set:construct_hash_table(test_set) for test_set in test_list}
        return test_structure
        
    def set_test_hash(self, current):
        '''
        current: current test list (gene list)
        '''        
        self.test_list = current
        return self.read_test_list()
    
    def read_test_list(self):
        '''
        load gene sets into memory

        differentiate between annotated and
        unannotated genes
        '''
        test_list = self.test_list
        if not exists(test_list):
            raise FileNotFoundError(
                ENOENT, strerror(ENOENT), test_list)
        else:
            t = loadtxt(test_list, dtype=str)
            self.unannotated_genes.append(set(t).difference(self.database_genes))
            return t
        
    def set_go_table(self):
        self.download_go()
        self.process_go()
        self.build()
        return
    
    def download_go(self):
        '''
        download GO annotations to file in current directory
        '''
        f = basename(self.link)
        # escape for Windows: \\
        go_file_path = '\\'.join([getcwd(), f])
        path_exists = exists(go_file_path)
        self.go_db = go_file_path
        
        if path_exists:
            pass
        else:
            print('Downloading GO annotation from file.')
            link = self.link
            link = get(link, stream = True)
            # write GO database to file in
            # memory conservative sizes
            with open(go_file_path, 'wb') as g:
                [g.write(chunk) for chunk in link.iter_content(chunk_size=1024) if chunk]
        return
        
    def process_go(self):
        '''
        convert database to data structure to use in testing
        '''
        go_db = self.go_db
        with open(go_db) as g: 
            file = [line.split('\t') for line in g]
        self.tmp_go_list = file
        return
        
    def build(self):
        '''
        Insert GO terms and associated genes into hash table
        '''
        print('Building GO term database.')
        [self.insert(_[7], _[5], _[0]) for _ in self.tmp_go_list if len(_) == 15]
        return
    
    def insert(self, aspect, go_term, iden):
        '''
        build hash table
        
        aspect: main Gene Ontology aspect
        go_term: GO term
        iden: gene identifier
        '''
        # set current parameters
        self.current_parameters = Param(aspect, go_term, iden)

        # insert identifer into global identifier set
        self.database_genes.add(iden)
        
        missing_aspect = self.missing_aspect()
        
        if missing_aspect:
            self.set_aspect()

        missing_go_term = self.missing_go_term()
        
        if missing_go_term:
            self.set_go_term()
            
        missing_iden = self.missing_iden()
       
        if missing_iden:
            self.set_iden()

        return
    
    def missing_aspect(self):
        '''
        decide whether to add aspect
        '''
        aspect = self.current_parameters.aspect
        return aspect not in self.go_table
    
    def missing_go_term(self):
        '''
        decide whether to add GO term
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        return go_term not in self.go_table[aspect]
    
    def missing_iden(self):
        '''
        decide whether to add identifer
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        iden = self.current_parameters.iden
        return iden not in self.go_table[aspect][go_term]       
    
    def set_aspect(self):
        '''
        set key to the newly- found aspect
        '''
        aspect = self.current_parameters.aspect
        self.go_table[aspect] = dict(population=list())
        return

    def set_go_term(self):
        '''
        set inner- most key to newly- found GO term
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        self.go_table[aspect][go_term] = list()
        return 
    
    def set_iden(self):
        '''
        set appended list element to newly- found GO identifier

        add identifier to total gene set
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        iden = self.current_parameters.iden
        self.go_table[aspect][go_term].append(iden)
        # add identifier to aspect population for statistical testing
        if iden not in self.go_table[aspect]['population']:
            self.go_table[aspect]['population'].append(iden)
        return

    def export_unannotated_genes(self):
        '''
        write unannotated genes to file
        '''
        for index, test_list in enumerate(self.test_list.keys()):
            output = test_list.replace('.tsv', '_unannotated.tsv')
            with open(output, 'w') as u:
                uw = writer(u, delimiter = '\n')
                uw.writerow(list(self.unannotated_genes[index]))
        return
