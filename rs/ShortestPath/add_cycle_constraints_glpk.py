from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Adds constraints to LP problem to eliminate cycles'

from copy import deepcopy

class AddCycleConstraintsGLPK(object):
    """
    Adds cycle constraints to ILP problem to prevent a shortest path with an identified cycle from
    being identified as a possible route to production for target compound, works with pyglpk
    package
    """
    def add_cycle_constraints(self, lp, allrxnsrev, cycle_constraints):
        '''
        Adds constraints to ILP to prevent a cycle from being the shortest path
        '''
        self.storecycleindexes = []
        self.storerxnindexes = {}
        self.totalrxnindexes = []
        self.lp = lp
        cyclestoic = []
        row_num = deepcopy(len(self.lp.rows))
        col_num = deepcopy(len(self.lp.cols))

        '''
        Adds number of variables needed for cycle constraints
        '''
        org_len = len(cyclestoic)
        temp = [0]*cycle_constraints.totalarcs
        cyclestoic = cyclestoic+temp
        new_len = len(cyclestoic)
        for i in range(org_len, new_len):
            self.storecycleindexes.append(i)

        self.lp.cols.add(len(cyclestoic))
        for i in range(col_num, len(self.lp.cols)):
            col = self.lp.cols[i]
            col.bounds = 0, 1
            col.name = 'Cycle Variable'
            col.kind = int

        for c in range(row_num, row_num+len(cyclestoic)):
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
        row_num = row_num-1
        for rxn in cycle_constraints.totalcyclerxns:
            for index in self.storerxnindexes[rxn]:
                row_num += 1
                self.lp.rows.add(1)
                r = self.lp.rows[row_num]
                rxnindex = allrxnsrev.index(rxn)
                r.matrix = [(rxnindex, -1), (col_num+index, 1)]
                r.bounds = 0, 0
                r.name = 'Cycle constraint'

        temp = []
        for i in self.storecycleindexes:
            temp.append((col_num+i, 1))
        row_num += 1
        self.lp.rows.add(1)
        r = self.lp.rows[row_num]
        r.matrix = temp
        r.bounds = 0, cycle_constraints.totalarcs-1
        r.name = 'Cycle constraint'
