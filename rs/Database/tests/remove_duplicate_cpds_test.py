from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that remove duplicate cpds in database'

import sqlite3
import sys
import glob
import os
import re
import unittest
import cobra
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import query as Q
from Database import remove_duplicate_cpds as rdc
PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)
init_db.Createdb(PATH+'/testinchi.db', True)
bkdb.BuildKbase(PATH+'/data7', PATH+'/KbasetoKEGGCPD_test.txt', PPATH+'/KbasetoKEGGRXN.txt', True, PATH+'/testinchi.db', 'bio')
DB = Q.Connector(PATH+'/testinchi.db')
DB.cnx.execute("UPDATE compound SET kegg_id=? WHERE ID=?", ('cpd00003', 'InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0'))
DB.conn.commit()

class RemoveDuplicateCpdsTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_removeduplicatecpds(self):
        rdc.OverlappingCpdIDs(PATH+'/testinchi.db')
        self.DB = Q.Connector(PATH+'/testinchi.db')
        S = self.DB.cnx.execute('select ID from compound')
        hits = S.fetchall()
        results = [i[0] for i in hits]
        self.assertFalse('cpd00003_c0' in results)

        S = self.DB.cnx.execute("select cpd_ID from reaction_compound where reaction_ID='%s'"%'rxn8_c0')
        hits = S.fetchall()
        results = [i[0] for i in hits]
        self.assertFalse('cpd00003_c0' in results)
        self.assertTrue('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0' in results)

        S = self.DB.cnx.execute("select cpd_ID from model_compound where model_ID='%s'" % 't2')
        hits = S.fetchall()
        results = [i[0] for i in hits]
        self.assertFalse('cpd00003_c0' in results)
        self.assertTrue('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0' in results)
        os.remove(PATH+'/testinchi.db')
if __name__ == '__main__':
    unittest.main() 

