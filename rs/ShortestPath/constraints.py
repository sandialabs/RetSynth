from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'retrieves necessary matricies for interger program'


import pulp
from copy import deepcopy
from tqdm import tqdm
import time
import re
from sys import platform

_OPTIMAL = True


def load_stoichometry_for_met(reactantrxns, productsrxns, allrxnsrev,
                              allrxnsrev_index, allrxnsrevset):
    '''Gets stoichometry for all compounds to be in the A matrix'''
    temp = {}
    react_not_in_all = reactantrxns.difference(allrxnsrevset)
    for s in react_not_in_all:
        if s + '_R' in allrxnsrev_index:
            temp[allrxnsrev_index[s + '_R']] = 1
        if s + '_F' in allrxnsrev_index:
            temp[allrxnsrev_index[s + '_F']] = -1
    react_in_all = reactantrxns.intersection(allrxnsrevset)
    product_not_in_all = productsrxns.difference(allrxnsrevset)
    for s in product_not_in_all:
        if s + '_R' in allrxnsrev_index:
            temp[allrxnsrev_index[s + '_R']] = -1
        if s + '_F' in allrxnsrev_index:
            temp[allrxnsrev_index[s + '_F']] = 1
    product_in_all = productsrxns.intersection(allrxnsrevset)
    for r in react_in_all:
        index = allrxnsrev_index[r]
        temp[index] = -1
    for p in product_in_all:
        index = allrxnsrev_index[p]
        temp[index] = 1
    return temp

class ConstructInitialLP(object):
    """Constructs A matrix and indidvidual reaction constraints"""
    def __init__(self, allrxns, allcpds, db, ignorerxns, gdbc, reverse_constraints=False):
        '''Initalize class'''
        self.allrxns = allrxns
        self.allcpds = deepcopy(allcpds)
        self.gdbc = gdbc
        self.DB = db
        self.reverse_constraints = reverse_constraints
        self.ignorerxns = ignorerxns
        self.A = []
        self.allcpds_new = []
        self.allrxnsrev = []
        self.allrxnsrevset = set()
        self.allrxnsrev_names = set()
        self.allrxnsrev_index = {}
        self.variables = []
        self.rxnnames = []
        self.lp = pulp.LpProblem('ShortestPath', pulp.LpMinimize)
        self.initial_reaction_constraints()
        if self.gdbc is True:
            self.allcpds = list(set(self.allcpds)) #ensures no compounds have duplicate entries
            self.initial_A_matrix(pulp)
        else:
            self.existing_A_matrix_pulp(pulp)

    def retrieve_stoichiometry(self, met):
        '''Retrieve stoichometry for each compound'''
        if met in self.reactant_dict:
            reactantrxns = self.reactant_dict[met]
        else:
            reactantrxns = set(self.DB.get_reactants_reactions(met))
            self.reactant_dict[met] = reactantrxns
        if met in self.product_dict:
            productsrxns = self.product_dict[met]
        else:
            productsrxns = set(self.DB.get_products_reactions(met))
            self.product_dict[met] = productsrxns
        stoich = load_stoichometry_for_met(reactantrxns, productsrxns,
                                           self.allrxnsrev,
                                           self.allrxnsrev_index,
                                           self.allrxnsrevset)
        return stoich

    def reaction_constraints_pulp(self, variable_name, rxn_name, pulp):
        '''Set reaction constraints (pulp)'''
        if rxn_name.startswith('EX_') or rxn_name in self.ignorerxns:
            variable = pulp.LpVariable(variable_name, cat=pulp.LpInteger,
                                       lowBound=0, upBound=0)
        else:
            variable = pulp.LpVariable(variable_name, cat=pulp.LpInteger,
                                       lowBound=0, upBound=1)
        self.variables.append(variable)

    def initial_reaction_constraints(self):
        '''Sets up column (individual reaction) constraints'''
        count = 0
        print ('STATUS: Generating reaction constraints ...')
        if self.reverse_constraints is False:
            for rxn in tqdm(self.allrxns):
                if self.DB.is_reversible_all(rxn) == 'true':
                    count += 1
                    self.rxnnames.append('R' + str(count))
                    self.allrxnsrev.append(str(rxn) + '_F')
                    self.allrxnsrevset.add(str(rxn) + '_F')
                    self.reaction_constraints_pulp('R' + str(count), str(rxn), pulp)

                    count += 1
                    self.rxnnames.append('R' + str(count))
                    self.allrxnsrev.append(str(rxn) + '_R')
                    self.allrxnsrevset.add(str(rxn) + '_R')
                    self.reaction_constraints_pulp('R' + str(count), str(rxn), pulp)

                else:
                    count += 1
                    self.rxnnames.append('R' + str(count))
                    self.allrxnsrev.append(str(rxn))
                    self.allrxnsrevset.add(str(rxn))
                    self.reaction_constraints_pulp('R' + str(count), str(rxn), pulp)

        else:
            for count, rxn in enumerate(self.allrxns):
                self.rxnnames.append('R' + str(count))
                self.allrxnsrev.append(str(rxn))
                self.allrxnsrevset.add(str(rxn))
                self.reaction_constraints_pulp(rxn + self.reverse_constraints[count], str(rxn), pulp)

        self.allrxnsrev_index = {key: index for index, key in enumerate(self.allrxnsrev)}

    def initial_A_matrix(self, solver=False):
        ''' Generates an matrix of compound constraints'''
        print('STATUS: Generating A matrix...')
        self.reactant_dict = {}
        self.product_dict = {}
        for count, met in enumerate(tqdm(self.allcpds)):
            stoich = self.retrieve_stoichiometry(met)
            if len(stoich) > 0:
                self.A.append(stoich)
                self.allcpds_new.append(met)
            else:
                pass
        self.allcpds = self.allcpds_new
        self.load_pulp_row_constraints(solver)

    def load_pulp_row_constraints(self, pulp):
        '''Loads constraints in to pulp integer linear problem'''
        print ('STATUS: Loading database A matrix (pulp)... ')
        for count, stoich in enumerate(tqdm(self.A)):
            self.lp += pulp.LpConstraint(pulp.lpSum(stoich[j]*self.variables[j] for j in stoich.keys()), name='c' + str(count) + ' constraint', sense=1, rhs=0)

    def existing_A_matrix_pulp(self, pulp):
        '''Loads existing matrix into pulp integer linear problem'''
        print ('STATUS: Generating compound constraints from preloaded file (pulp) ...')
        for count, stoich in enumerate(tqdm(self.gdbc)):
            self.lp += pulp.LpConstraint(pulp.lpSum(stoich[j]*self.variables[j] for j in stoich.keys()), name='c' + str(count) + ' constraint', sense=1, rhs=0)
