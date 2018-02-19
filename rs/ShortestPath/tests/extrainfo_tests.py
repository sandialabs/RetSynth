from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests extrainfo_tests.py'

import unittest
import re
import os
from ShortestPath import constraints as co
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from ShortestPath import extractinfo as ei
PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

'''CONNECT TEST DATABASE'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data2', '../../Database/KbasetoKEGGCPD.txt',
                '../../Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
inmets = DB.get_compounds_in_model('t1')
osp=[['rxn1_c0_F', 'rxn5_c0', 'rxn6_c0', 'rxn7_c0'],
     ['rxn1_c0_F', 'rxn2_c0', 'rxn3_c0', 'rxn4_c0']]
os.remove(PATH+'/test.db')

class ExtrainfoTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_extrainfo(self):
        """Tests A matrix construction glpk"""
        print ("Tests that correct information about compounds,\
                       reactions and organisms is being retrieved")

        ex_info = ei.Extract_Information(osp, inmets, DB)
        os_dict = ex_info.temp_rxns[1]
        external_metabolites = ex_info.temp_exmets[1]

        print ("...Testing that correct number of external compounds was retrieved")
        self.assertEqual(len(external_metabolites), 5)
        self.assertEqual(len(os_dict.keys()), 4)
 
        print ("...Testing that correct reaction information (directionality and name)\
                       was retrieved")
        self.assertEqual(len(os_dict['rxn1_c0'].keys()), 6)
        self.assertEqual(os_dict['rxn1_c0']['direction'], 'forward')
        self.assertEqual(os_dict['rxn1_c0']['name'], 'rxn1_c0')

        print ("...Testing that correct reactants were retrieved")
        self.assertEqual(len(os_dict['rxn1_c0']['reactants']), 1)
        self.assertIn('cpdA_c0', os_dict['rxn1_c0']['reactants'])

        print ("...Testing that correct products were retrieved")
        self.assertEqual(len(os_dict['rxn1_c0']['products']), 2)
        self.assertIn('cpdT_c0', os_dict['rxn1_c0']['products'])
        self.assertIn('cpdBy1_c0', os_dict['rxn1_c0']['products'])

if __name__ == '__main__':
    unittest.main()
