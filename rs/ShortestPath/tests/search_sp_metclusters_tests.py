from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests search paths through all metabolic clusters'
import re
import os
import unittest
from ShortestPath import constraints as co
from Database import query as Q
from ShortestPath import search_sp_metclusters as smc
from ShortestPath import integerprogram_glpk as ip_glpk
from ShortestPath import integerprogram_pulp as ip_pulp
from Parser import generate_output as go
from Database import generate_database as gen_db

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/ShortestPath/tests', '', PATH)

'''CONNECT TEST DATABASE'''
gen_db.Createdb(PATH+'/test.db', PATH+'/data3', False, 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
OUTPUT = go.Output(DB, PATH, False, False)
os.remove(PATH+'/test.db')

class ExtrainfoTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_search_sp_metclusters_glpk(self):
        """Tests search for pathway for a given compound in all metabolic clusters, glpk"""
        print ('...Tests search for pathway for a given compound in all metabolic clusters, glpk')
        try:
            import glpk
            LP = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'GLPK')
            IP = ip_glpk.IntergerProgram(DB, 10, 20, 0, True)
            S = smc.SearchMetabolicClusters('cpdT_c0', LP, IP, OUTPUT, DB)
            self.assertEqual(len(S.total_sp_clusters), 2)
            self.assertIn(['cpdT_c0', 4, ['t1', 't3']], S.total_sp_clusters)
            if os.path.isfile(PPATH+'/path_length_all_organism_cpdT_c0.txt') is True:
                os.remove(PPATH+'/path_length_all_organism_cpdT_c0.txt')
            if os.path.isfile(PATH+'/path_length_all_organism_cpdT_c0.txt') is True:
                os.remove(PATH+'/path_length_all_organism_cpdT_c0.txt')
            if os.path.isfile(PPATH+'/path_length_all_organism_None.txt') is True:
                os.remove(PPATH+'/path_length_all_organism_None.txt')
            if os.path.isfile(PATH+'/path_length_all_organism_None.txt') is True:
                os.remove(PATH+'/path_length_all_organism_None.txt')                
            if os.path.isfile(PPATH+'/optimal_pathways.txt') is True:
                os.remove(PPATH+'/optimal_pathways.txt')
            if os.path.isfile(PATH+'/optimal_pathways.txt') is True:
                os.remove(PATH+'/optimal_pathways.txt')
        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_search_sp_metclusters_pulp(self):
        """Tests search for pathway for a given compound in all metabolic clusters, pulp"""
        print ('...Tests search for pathway for a given compound in all metabolic clusters, pulp')
        try:
            import pulp
            LP = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'PULP')
            IP = ip_pulp.IntergerProgram(DB, 10, 20, 0,True)
            S = smc.SearchMetabolicClusters('cpdT_c0', LP, IP, OUTPUT, DB)
            self.assertEqual(len(S.total_sp_clusters), 2)
            self.assertIn(['cpdT_c0', 4, ['t1', 't3']], S.total_sp_clusters)
            if os.path.isfile(PPATH+'/path_length_all_organism_cpdT_c0.txt') is True:
                os.remove(PPATH+'/path_length_all_organism_cpdT_c0.txt')
            if os.path.isfile(PATH+'/path_length_all_organism_cpdT_c0.txt') is True:
                os.remove(PATH+'/path_length_all_organism_cpdT_c0.txt')
            if os.path.isfile(PPATH+'/path_length_all_organism_None.txt') is True:
                os.remove(PPATH+'/path_length_all_organism_None.txt')
            if os.path.isfile(PATH+'/path_length_all_organism_None.txt') is True:
                os.remove(PATH+'/path_length_all_organism_None.txt')
            if os.path.isfile(PPATH+'/optimal_pathways.txt') is True:
                os.remove(PPATH+'/optimal_pathways.txt')
            if os.path.isfile(PATH+'/optimal_pathways.txt') is True:
                os.remove(PATH+'/optimal_pathways.txt')
        except ImportError:
            print ('pulp package not installed, cannot run tests for this python package')

if __name__ == '__main__':
    unittest.main()
