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
    def __init__(self, allrxns, allcpds, db, ignorerxns, includerxns, forward_direction=True, reverse_direction=False,
                 lp=None, variables=None, allrxnsrev_dict_rev=None,
                 allrxnsrev_dict=None, allrxnsrev=None):
        '''Initalize class'''

        self.allrxns = allrxns
        self.allcpds = deepcopy(allcpds)
        self.DB = db
        self.ignorerxns = ignorerxns
        self.includerxns = includerxns
        self.forward_direction = forward_direction
        self.reverse_direction = reverse_direction
        # self.lp = pulp.LpProblem('ShortestPath', pulp.LpMinimize)

        if lp is  None:
            self.A = []
            self.allcpds_new = []
            self.allrxnsrev = []
            self.allrxnsrev_dict = {}
            self.allrxnsrev_dict_rev = {}
            self.allrxnsrev_names = set()
            self.allrxnsrev_index = {}
            self.variables = []
            self.variables_load = []
            self.lp = pulp.LpProblem('ShortestPath', pulp.LpMinimize)
            self.allcpds = list(set(self.allcpds)) #ensures no compounds have duplicate entries
            
            self.initial_reaction_constraints()
            self.initial_A_matrix(pulp)
           
            if self.ignorerxns:
                self.reaction_constraints_ignore_reactions()

            if self.includerxns:
                self.reaction_constraints_include_reactions()


        else:
            print ('STATUS: Building pre-stored variables...')
            self.lp = lp
            self.variables = self.lp.variables()
            self.allrxnsrev_dict_rev = allrxnsrev_dict_rev
            self.allrxnsrev_dict = allrxnsrev_dict
            self.allrxnsrev = allrxnsrev

            if self.ignorerxns:
                self.reaction_constraints_ignore_reactions()

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
        self.variables = self.lp.variables()

    ###GENERATE REACTION VARIABLE CONSTRAINTS
    def reaction_constraints_pulp(self, variable_name, rxn_name, pulp):
        '''Set reaction constraints (pulp)'''
        if rxn_name.startswith('EX_'):
            variable = pulp.LpVariable(variable_name, cat=pulp.LpInteger,
                                       lowBound=0, upBound=0)
        else:
            variable = pulp.LpVariable(variable_name, cat=pulp.LpInteger,
                                       lowBound=0, upBound=1)
        self.variables_load.append(variable)

    def _update_reaction_constraints_for_ingorerxns(self, reaction_id):
        tmp = [k for k in self.variables if k.name==reaction_id]
        variable = tmp[0]
        variable.upBound = 0 
        variable.lowBound = 0

    def _update_reaction_constraints_for_includerxns(self, reaction_id):
        tmp = [k for k in self.variables if k.name==reaction_id]
        variable = tmp[0]
        variable.upBound = 1 
        variable.lowBound = 1
 
 
    def reaction_constraints_ignore_reactions(self):
        for reaction in self.ignorerxns:
            try:
                reaction_id = self.allrxnsrev_dict_rev[reaction]
                self._update_reaction_constraints_for_ingorerxns(reaction_id) 
            except KeyError:
                try:
                    reaction_id_forward = self.allrxnsrev_dict_rev[reaction+'_F']
                    reaction_id_reverse = self.allrxnsrev_dict_rev[reaction+'_R']
                    self._update_reaction_constraints_for_ingorerxns(reaction_id_forward) 
                    self._update_reaction_constraints_for_ingorerxns(reaction_id_reverse) 

                except:
                    print ('WARNING: {} is not in reaction list'.format(reaction))


    def reaction_constraints_include_reactions(self):
        for reaction in self.includerxns:
            try:
                reaction_id = self.allrxnsrev_dict_rev[reaction]
                self._update_reaction_constraints_for_ingorerxns(reaction_id) 
            except KeyError:
                try:
                    reaction_id_forward = self.allrxnsrev_dict_rev[reaction+'_F']
                    reaction_id_reverse = self.allrxnsrev_dict_rev[reaction+'_R']

                    if self.forward_direction:
                        self._update_reaction_constraints_for_includerxns(reaction_id_forward)

                    elif self.reverse_direction:
                        self._update_reaction_constraints_for_includerxns(reaction_id_reverse)
                    
                    elif self.reverse_direction and self.forward_direction:
                        print ('WARNING: cannot have both directions of a reaction be apart of a pathway as this would create a cycle...continuing to solve without forcing the specified reaction to be apart of pathway')

                except:
                    print ('WARNING: {} is not in reaction list'.format(reaction))


    def load_reaction_variables(self, rxn_rev, rxn, rxn_id):
        self.allrxnsrev.append(rxn_rev)
        self.reaction_constraints_pulp(rxn_id, str(rxn), pulp)
        self.allrxnsrev_dict[rxn_id] = rxn_rev
        self.allrxnsrev_dict_rev[rxn_rev] = rxn_id

    def initial_reaction_constraints(self):
        '''Sets up column (individual reaction) constraints'''
        count = 0
        print ('STATUS: Generating reaction constraints ...')
        for rxn in tqdm(self.allrxns):
            if self.DB.is_reversible_all(rxn) == 'true':
                count += 1
                self.load_reaction_variables(str(rxn) + '_F', str(rxn), 'R' + str(count))

                count += 1
                self.load_reaction_variables(str(rxn) + '_R', str(rxn), 'R' + str(count))

            else:
                count += 1
                self.load_reaction_variables(str(rxn), str(rxn), 'R' + str(count))

        self.allrxnsrevset = set(self.allrxnsrev)
        self.allrxnsrev_index = {key: index for index, key in enumerate(self.allrxnsrev)}
        self.variables = deepcopy(self.variables_load)