from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests cycle.py and addition of cycle constraints,\
                   ensures it obtains correct answers'

import unittest
import re
import os
from copy import deepcopy
from ShortestPath import constraints as co
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from ShortestPath import integerprogram_glpk as ip_glpk
from ShortestPath import integerprogram_pulp as ip_pulp
from ShortestPath import cycle as cy
from ShortestPath import add_cycle_constraints_glpk as acc_glpk
from ShortestPath import add_cycle_constraints_pulp as acc_pulp

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/ShortestPath/tests', '', PATH)

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

'''GENERATE TEST DATABASE'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data2', PPATH+'/Database/KbasetoKEGGCPD.txt',
                PPATH+'/Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
os.remove(PATH+'/test.db')

class CycleTests(unittest.TestCase):
    """Test class"""
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.CYCLE = cy.CycleCheck(DB)

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
        self.CYCLE

    def test_identification_of_cycle_glpk(self):
        """Tests identification of cycles with glpk package"""
        try:
            import glpk
            print ("...Testing identification of cycles with glpk package")
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'GLPK')
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 0, 'False')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            for solution in osp:
                cycletest = self.CYCLE.run_cycle_check(solution, inmets)
                self.assertEqual(cycletest, True)
                print ("...Testing that cycle variables and arcs are correct")
                self.assertEqual(len(self.CYCLE.totalvariables), 9)
                self.assertEqual(self.CYCLE.totalarcs, 9)

        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_identification_of_k_path_glpk(self):
        """Tests identification of cycles with glpk package"""
        try:
            import glpk
            print ("...Testing identification of cycles with glpk package")
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'GLPK')
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 1, 'False')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            self.assertEqual(len(osp), 3)
            self.assertEqual(len(osp[-1]), 6)

        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_identification_of_cycle_pulp(self):
        """Tests identification of cycles with pulp package"""
        try:
            import pulp
            print ("...Testing identification of cycles with pulp package")
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'PULP')
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 0, 'False', 'None')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            for solution in osp:
                cycletest = self.CYCLE.run_cycle_check(solution, inmets)
                self.assertEqual(cycletest, True)
                print ("...Testing that cycle variables and arcs are correct")
                self.assertEqual(len(self.CYCLE.totalvariables), 9)
                self.assertEqual(self.CYCLE.totalarcs, 9)

        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')

    def test_identification_of_k_path_pulp(self):
        """Tests identification of cycles with pulp package"""
        try:
            import pulp
            print ("...Testing identification of cycles with pulp package")
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'PULP')
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 1, 'False', 'None')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            self.assertEqual(len(osp), 3)
            self.assertEqual(len(osp[-1]), 6)

        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')

    def test_addition_of_cycle_constraints_glpk(self):
        """Tests the addition of cycle constraints for glpk package"""
        try:
            import glpk
            print ("...Testing the addition of cycle constraints for glpk package")
            ACC = acc_glpk.AddCycleConstraintsGLPK()
            print ("Tests if cycle constraints are getting added correctly")
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            LP = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'GLPK')
            original_row = deepcopy(len(LP.lp.rows))
            original_col = deepcopy(len(LP.lp.cols))
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 0, 'False')
            osp = IP.run_glpk(LP, inmets, inrxns, 'cpdT_c0', 'True')
            count_arcs = 0
            for solution in osp:
                cycletest = self.CYCLE.run_cycle_check(solution, inmets)
                ACC.add_cycle_constraints(LP.lp, LP.allrxnsrev, self.CYCLE)
                self.assertEqual(self.CYCLE.totalarcs, len(ACC.storecycleindexes))
                count_arcs += self.CYCLE.totalarcs
            self.assertEqual(original_col+count_arcs, len(ACC.lp.cols))
            self.assertEqual(original_row+count_arcs+len(osp), len(ACC.lp.rows))
        except ImportError:
            print ('glpk package not installed, cannot run tests for this python package')

    def test_addition_of_cycle_constraints_pulp(self):
        """Tests the addition of cycle constraints for pulp package"""
        try:
            import pulp
            print ("...Testing the addition of cycle constraints for pulp package")
            ACC = acc_pulp.AddCycleConstraints()
            print ("Tests if cycle constraints are getting added correctly")
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            LP = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'PULP')
            original_row = deepcopy(len(LP.lp.constraints))
            original_col = deepcopy(len(LP.variables))
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 0, 'False', 'None')
            osp = IP.run_glpk(LP, inmets, inrxns, 'cpdT_c0', 'True')
            count_arcs = 0
            for solution in osp:
                cycletest = self.CYCLE.run_cycle_check(solution, inmets)
                ACC.add_cycle_constraints(LP.lp, LP.variables, LP.allrxnsrev, self.CYCLE)
                LP.variables = ACC.variables
                self.assertEqual(self.CYCLE.totalarcs, len(ACC.storecycleindexes))
                count_arcs += self.CYCLE.totalarcs
            self.assertEqual(original_col+count_arcs, len(ACC.variables))
            self.assertEqual(original_row+count_arcs+len(osp), len(ACC.lp.constraints))
        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')
if __name__ == '__main__':
    unittest.main()
