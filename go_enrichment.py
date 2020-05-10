from set_parameters import Param
from copy import deepcopy
from scipy.stats import hypergeom
from re import sub
from csv import writer

class G():
    def __init__(self, t=None, s=None):
        '''
        t: Struct class instance
        s: significance threshold
        '''
        self.test_list = t.test_list
        self.go_dict = t.go_dict
        self.unannotated_genes = t.unannotated_genes
        self.significance = s
        self.hypergeom_test()
        self.export_results()
        
    def hypergeom_test(self):
        '''
        set hyper- geometric testing parameters and conduct test
        '''
        if not self.go_dict:
            raise Warning('Missing Gene Ontology.')
        
        self.struct = deepcopy(self.go_dict)
        test_list = self.test_list

        if not test_list:
            raise Warning('Missing testing list.')

        test_set_hash_table = self.test_dict_import
        print('Hypergeometric testing.')
        self.struct = {name:test_set_hash_table(test_set) for name,test_set in test_list.items()}
        return 
            
    def test_dict_import(self, current):
        '''
        add unannotated genes to unannotated set

        import test set and conduct test       
        
        current: gene list to analyze
        '''
        # set current test list
        self.test_list = current
        return self.collect_parameters_and_conduct_test()
    
    def collect_parameters_and_conduct_test(self):
        '''
        test for over- representation between test set and database

        filter results on significant value provided
        '''
        for aspect,go_term_set in self.struct.items():
            self.set_population_parameters()
            for go_term,identifiers in go_term_set.items():
                
                # set current parameters
                self.current_parameters = Param(aspect, go_term, identifiers)
                
                # check if accessing the aspects's population parameters
                self.population_element = (self.current_parameters.go_term == 'population')
                self.edit_struct()
                if self.population_element:
                    pass
                else:
                    # test for over- representated with hyper-
                    # geometric test
                    self.over_representation_test()

                    # identify terms with significant enrichment
                    iden = self.current_parameters.iden
                    current = self.struct[aspect][go_term]
                    test_value = current[len(current)-1]
                    if test_value <= self.significance:
                        self.struct[aspect][go_term] = iden
                        self.struct[aspect][go_term].insert(0, test_value)

                    else:
                        self.struct[aspect][go_term] = [test_value, 'NA']
        return self.struct
            
    def set_population_parameters(self):
        '''
        population parameter constructor
        '''
        self.N, self.k = None,None
    
    def edit_struct(self):
        '''
        over- write GO database copy with hyper- geometric
        testing parameters
        '''
        self.set_vanilla_struct()
        self.set_overlap_size()
        self.set_size()
        if self.population_element:
            pass
        else:
            # append population parameters to
            # sample parameter list
            self.simplify_struct()
        return
    
    def set_vanilla_struct(self):
        '''
        reset identifier list to an empty list
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        self.struct[aspect][go_term] = list()
        return
    
    def set_overlap_size(self):
        '''
        get number of successes in population
        and sample
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        k = self.intersection()
        if self.population_element:
            self.k = k
        self.struct[aspect][go_term].append(k)
        return
    
    def intersection(self):
        '''
        return intersection size for current test
        list and identifier list

        identify genes in test list without an
        annotation
        '''
        iden_list = set(self.current_parameters.iden)
        test_list = set(self.test_list)
        intersect = list(iden_list & test_list)
        return len(intersect)
    
    def set_size(self):
        '''
        get number of successes and failures
        in population and sample
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        size = len(self.go_dict[aspect][go_term])
        if self.population_element:
            self.N = size
        self.struct[aspect][go_term].append(size)
        return 
        
    def simplify_struct(self):
        '''
        insert testing parameters for identifer
        list to apply hyper- geometric test
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        N, k = self.N, self.k
        self.struct[aspect][go_term].insert(1, N)
        self.struct[aspect][go_term].insert(2, k)
        return
    
    def over_representation_test(self):
        '''
        apply hyper- geometric test

        only append value if less than significance
        '''
        aspect = self.current_parameters.aspect
        go_term = self.current_parameters.go_term
        current = self.struct[aspect][go_term]
        result = [hypergeom.sf(_[0]-1, _[1], _[2], _[3]) for _ in [current]]
        result = result.pop()
        current.append(result)
        return
    
    def export_results(self):
        '''
        export aspect, go_term and significance: n^3 complexity

        export unannotated genes from test list(s)
        '''
        for test_list,enrichment in self.struct.items():
            output = test_list.replace('.tsv', '_enrichment.csv')
            with open(output, 'a') as current_file:                
                current_file.write("%s,%s,%s,%s\n"%('aspect', 'go_term', 'significance', 'identifiers'))
                for aspect, go_terms in enrichment.items():
                    for terms, iden in go_terms.items():
                        if terms != 'population':
                            current_file.write("%s,%s,%s,%s\n"%(aspect, terms, iden[0], ' '.join(iden[1:])))

        with open('unannotated.csv', 'w') as f:
            u = writer(f, delimiter = '\n')
            u.writerow(list(self.unannotated_genes[0])                
        return


    
