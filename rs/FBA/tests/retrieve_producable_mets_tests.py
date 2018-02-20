from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests optimizing of a FBA model from the Database \
                   and addition of external reactions to model'
import os
import re
import unittest
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from FBA import retrieve_producable_mets as rpm

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/FBA/tests', '', PATH)

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

'''Load database'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', PPATH+'/Database/KbasetoKEGGCPD.txt',
                PPATH+'/Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
inmets = DB.get_compounds_in_model('t3')
inrxns = DB.get_reactions_in_model('t3')
os.remove(PATH+'/test.db')

class Retrieve_producable_compoundsTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.R = rpm.RetrieveActiveRxnsCompounds('t3', inmets, inrxns, DB, False)

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
        self.R
    def test_retrieve_producable_compounds(self):
        """tests to see if code can accurately identify
           compounds that can be produced with media"""
        print ("""Tests to see if code can accurately identify \
                      compounds that can be produced with media""")
        self.assertIn('cpdI_e0', self.R['t3'][0])
        self.assertIn('cpdI_c0', self.R['t3'][0])
        self.assertIn('cpdH_c0', self.R['t3'][0])
        self.assertIn('cpdG_c0', self.R['t3'][0])
        self.assertIn('cpdF_c0', self.R['t3'][0])
        self.assertFalse('cpdA_c0' in self.R['t3'][0])

if __name__ == '__main__':
    unittest.main()
