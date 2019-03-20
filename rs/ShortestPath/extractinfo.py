from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Pull out all reactants and products for reactions \
                  that need to be inserted to organism'

import re
from copy import deepcopy

class Extract_Information(object):
    """
    Retrieves names of compounds, reactions and organisms for
    reactions and compounds that need to be added
    to produce a target compound
    """
    def __init__(self, optimal_pathways, incpds, inrxns, db):
        '''Initalize class'''
        self.temp_rxns = {}
        self.temp_exmets = {}
        self.temp_external = {}
        self.inrxns = inrxns
        self.DB = db
        self.incpds = incpds
        for count, path in enumerate(optimal_pathways):
            count += 1
            if path:
                os_dict, excpds, count_external = self.extractinfo(path)
                if os_dict is not None:
                    self.temp_rxns[count] = os_dict
                    self.temp_exmets[count] = excpds
                    self.temp_external[count] = count_external
    def get_info(self, rxn, path_dict, excpds, Direction=False):
        '''
        Gets information for reactions, compounds and organisms
        '''
        path_dict[rxn]['organisms'] = []
        path_dict[rxn]['genes'] = []
        for org in self.DB.get_reaction_species(rxn):
            path_dict[rxn]['organisms'].append(org)
            path_dict[rxn]['genes'].append(self.DB.get_genes(rxn, org))
        list(set(path_dict[rxn]['genes']))

        path_dict[rxn]['name'] = self.DB.get_reaction_name(rxn)
        path_dict[rxn]['reactants'] = {}
        path_dict[rxn]['products'] = {}
        if Direction is True:
            for prod in self.DB.get_products(rxn):
                path_dict[rxn]['products'][prod] = self.DB.get_compound_name(prod)
                if prod not in self.incpds:
                    excpds.append(prod)
            for react in self.DB.get_reactants(rxn):
                path_dict[rxn]['reactants'][react] = self.DB.get_compound_name(react)
                if react not in self.incpds:
                    excpds.append(react)
        else:
            for react in self.DB.get_reactants(rxn):
                path_dict[rxn]['reactants'][react] = self.DB.get_compound_name(react)
                if react not in self.incpds:
                    excpds.append(react)
            for prod in self.DB.get_products(rxn):
                path_dict[rxn]['products'][prod] = self.DB.get_compound_name(prod)
                if prod not in self.incpds:
                    excpds.append(prod)
        return(path_dict, excpds)


    def extractinfo(self, path):
        '''
        Builds necessary dictionaries to hold the info found in the get_info function
        '''
        path_dict = {}
        excpds = []
        count = 0
        count_external = 0
        for rxn in path:
            rxn = deepcopy(rxn)
            if re.search('_R$', rxn) is not None:
                rxn = re.sub('_R$', '', rxn)
                if rxn not in self.inrxns:
                    count_external+=1
                if rxn not in path_dict.keys():
                    path_dict[rxn] = {}
                    path_dict[rxn]['direction'] = 'reverse'
                    path_dict, excpds = self.get_info(rxn, path_dict, excpds, Direction=True)
                else:
                    count += 1
            else:
                if re.search('_F$', rxn) is not None:
                    rxn = re.sub('_F$', '', rxn)
                if rxn not in self.inrxns:
                    count_external+=1
                if rxn not in path_dict.keys():
                    path_dict[rxn] = {}
                    path_dict[rxn]['direction'] = 'forward'
                    path_dict, excpds = self.get_info(rxn, path_dict, excpds, Direction=False)
                else:
                    count += 1
        excpds = list(set(excpds))
        if count == 0:
            return(path_dict, excpds, count_external)
        else:
            print ('WARNING: Solution may have reversible reactions')
            return(path_dict, excpds, count_external)
