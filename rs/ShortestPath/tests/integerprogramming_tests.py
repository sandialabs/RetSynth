from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests constraints.py (which builds the A matrix and reaction bounds) \
                   with a test database'

import unittest
import re
import os
from ShortestPath import constraints as co
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from ShortestPath import integerprogram_glpk as ip_glpk
from ShortestPath import integerprogram_pulp as ip_pulp
PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/ShortestPath/tests', '', PATH)

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', PPATH+'/Database/KbasetoKEGGCPD.txt',
                PPATH+'/Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
os.remove(PATH+'/test.db')


class ConstraintsTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_integerprograms_glpk(self):
        """Tests actual IP solver, glpk"""
        print ("Tests actual IP solver, glpk")
        try:
            import glpk
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'GLPK')
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 0, 'False')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'False')
            self.assertEqual(len(osp), 1)
            self.assertEqual(len(osp[0]), 4)
            self.assertIn('rxn1_c0_F', osp[0])
        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_integerprograms_pulp(self):
        """Tests actual IP solver, pulp"""
        print ("Tests actual IP solver, pulp")
        try:
            import pulp
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'PULP')
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 0, 'False', 'None')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'False')
            self.assertEqual(len(osp), 1)
            self.assertEqual(len(osp[0]), 4)
            self.assertIn('rxn1_c0_F', osp[0])
        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')

    def test_run_glpk_multipaths_glpk(self):
        """Tests if multiple pathways can be identified glpk"""
        print ("Tests if multiple pathways can be identified")
        try:
            import glpk
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True,
                                      reverse_constraints=False, specified_pysolver='GLPK')
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 0, 'False')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            self.assertEqual(len(osp), 2)
            self.assertEqual(len(osp[0]), 4)
        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_run_glpk_multipaths_pulp(self):
        """Tests if multiple pathways can be identified pulp"""
        print ("Tests if multiple pathways can be identified pulp")
        try:
            import pulp
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True,
                                      reverse_constraints=False, specified_pysolver='PULP')
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 0, 'False', 'None')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            self.assertEqual(len(osp), 2)
            self.assertEqual(len(osp[0]), 4)
        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')

    def test_run_glpk_multipaths_cycle_glpk(self):
        """Tests if multiple pathways can be identified glpk"""
        print ("Tests if multiple pathways can be identified")
        try:
            import glpk
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True,
                                      reverse_constraints=False, specified_pysolver='GLPK')
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 0, 'True')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            self.assertEqual(len(osp), 1)
            print (osp)
            self.assertEqual(len(osp[0]), 5)
        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_run_glpk_multipaths_cycle_pulp(self):
        """Tests if multiple pathways can be identified pulp"""
        print ("Tests if multiple pathways can be identified pulp")
        try:
            import pulp
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True,
                                      reverse_constraints=False, specified_pysolver='PULP')
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 0, 'True', 'None')
            inmets = DB.get_compounds_in_model('t1')
            inrxns = DB.get_reactions_in_model('t1')
            osp = IP.run_glpk(C, inmets, inrxns, 'cpdT_c0', 'True')
            self.assertEqual(len(osp), 1)
            self.assertEqual(len(osp[0]), 5)
        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')

if __name__ == '__main__':
    unittest.main()
