from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Adds constraints to LP problem to eliminate cycles'

import sys
import itertools
from copy import deepcopy
import pulp

class AddCycleConstraints(object):
    """
    Adds cycle constraints to ILP problem to prevent a shortest path with an identified cycle from
    being identified as a possible route to production for target compound, works with pyglpk package
    """
    def add_cycle_constraints(self, lp, variables, allrxnsrev, cycle_constraints):
        '''
        Adds constraints to ILP to prevent a cycle from being the shortest path
        '''

        '''Set initial variables'''
        self.storecycleindexes = []
        self.storerxnindexes = {}
        self.totalrxnindexes = []
        self.variables = variables
        self.lp = lp
        self.cycle_variables = []
        cyclestoic = []
        row_num = len(self.lp.constraints)

        '''
        Adds number of variables needed for cycle constraints
        '''
        org_len = len(cyclestoic)
        temp = [0]*cycle_constraints.totalarcs
        cyclestoic = cyclestoic+temp
        new_len = len(cyclestoic)
        for i in range(org_len, new_len):
            self.storecycleindexes.append(i)

        for c in range(len(self.variables), len(self.variables)+len(cyclestoic)):
            variable = pulp.LpVariable('cycle '+str(c), cat=pulp.LpInteger, lowBound=0, upBound=1)
            self.cycle_variables.append(variable)
            self.totalrxnindexes.append(c)

        '''
        Specifies which added variables correspond to rxn and compound arc
        '''
        count = 0
        for rxn in cycle_constraints.totalcyclerxns:
            self.storerxnindexes[rxn] = []
            for i in range(count, cycle_constraints.totalcyclerxns[rxn]+count):
                self.storerxnindexes[rxn].append(self.storecycleindexes[i])
            count += cycle_constraints.totalcyclerxns[rxn]

        '''
        Give values to the new constraints
        '''
        for rxn in cycle_constraints.totalcyclerxns:
            for index in self.storerxnindexes[rxn]:
                row_num += 1
                temp = [0]*len(cyclestoic)
                rxnindex = allrxnsrev.index(rxn)
                self.lp += pulp.LpConstraint((-1*variables[rxnindex]+1*self.cycle_variables[index]), name='Cycle constraint '+str(row_num), sense=0, rhs=0)

        self.variables = self.variables+self.cycle_variables
        row_num += 1
        self.lp += pulp.LpConstraint(pulp.lpSum(1*self.cycle_variables[i] for i in range(len(self.cycle_variables))), sense=-1, name='Cycle constraint '+str(row_num), rhs=cycle_constraints.totalarcs-1)
