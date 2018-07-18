from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Checks for cycles in answers'

import re

def verbose_print(verbose, line):
    if verbose:
        print(line)

class CycleCheck(object):
    """Runs check to see if a cycle is present in solution (shortest path)"""
    def __init__(self, db, verbose):
        '''Initalize class'''
        self.DB = db
        self.verbose = verbose

    def run_cycle_check(self, optimal_pathway, incpds):
        '''
        Run cycle check
        '''
        self.totalarcs = []
        self.totalvariables = []
        totalvariables = []
        totalarcs = 0
        self.totalcyclerxns = {}
        for rxn in optimal_pathway:
            totalvariables.append(rxn)
            match = re.search('_F$', rxn)
            rxn_org = rxn
            if match is not None:
                rxn_org = re.sub('_F$', '', rxn)
            match = re.search('_R$', rxn)
            if match is not None:
                rxn_org = re.sub('_R$', '', rxn)
            self.totalcyclerxns[rxn] = 0
            for reactant in self.DB.get_reactants(rxn_org):
                if reactant  not in incpds:
                    totalvariables.append(reactant)
                    totalarcs += 1
                    self.totalcyclerxns[rxn] += 1
            for product in self.DB.get_products(rxn_org):
                if product not in incpds:
                    totalvariables.append(product)
                    self.totalcyclerxns[rxn] += 1
                    totalarcs += 1
        if totalarcs > len(set(totalvariables))-1:
            self.totalvariables = list(set(totalvariables))
            self.totalarcs = totalarcs
        if self.totalvariables:
            verbose_print(self.verbose, 'STATUS: optimal pathway has cycle')
            return (True)
        else:
            verbose_print(self.verbose, 'STATUS: No cycles were found')
            return (False)
