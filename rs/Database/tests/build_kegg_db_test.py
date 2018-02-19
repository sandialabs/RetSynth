from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that add reactions and compounds to database from kegg'

import os
import re
import sqlite3
import unittest
from copy import deepcopy
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import build_KEGG_db as bkeggdb
PATH = os.path.dirname(os.path.abspath(__file__))
init_db.Createdb(PATH+'/kbasetestadd.db', False)
bkdb.BuildKbase(PATH+'/datam', '../KbasetoKEGGCPD.txt', '../KbasetoKEGGRXN.txt', False, PATH+'/kbasetestadd.db', 'bio')
DB = Q.Connector(PATH+'/kbasetestadd.db')
kbaserxnsadd = DB.get_all_reactions()
kbaserxnsadd1 = deepcopy(kbaserxnsadd)

init_db.Createdb(PATH+'/kbasetestaddinchi.db', True)
bkdb.BuildKbase(PATH+'/data', '../KbasetoKEGGCPD.txt', '../KbasetoKEGGRXN.txt', True, PATH+'/kbasetestaddinchi.db', 'bio')
DBinchi = Q.Connector(PATH+'/kbasetestaddinchi.db')
kbaserxnsinchi = DBinchi.get_all_reactions()
kbaserxns1inchi = deepcopy(kbaserxnsinchi)

class BuildKEGG(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_KEGG_with_only_KEGG_orgs(self):
        """test KEGG database construction"""
        K = bkeggdb.CompileKEGGIntoDB(PATH+'/kbasetestadd.db', 'bacteria', False, 4, 2, 2, 'bio', True)
        DBadd = Q.Connector(PATH+'/kbasetestadd.db')
        newrxns = set()
        with open('newrxns.txt') as fin:
            for line in fin:
                line = line.strip()
                newrxns.add(line)

        inter = set(newrxns) - set(kbaserxnsadd1)
        kbaserxnsadd = DBadd.get_all_reactions()
        self.assertEqual(len(kbaserxnsadd1)+len(inter), len(kbaserxnsadd))

        reactants = DBadd.get_reactants('R02187_c0')
        products = DBadd.get_products('R02187_c0')
        self.assertIn('C00404_c0', reactants)
        self.assertIn('C00221_c0', reactants)
        self.assertIn('C00404_c0', products)
        self.assertIn('C01172_c0', products)

        reactants = DBadd.get_reactants('R07618_c0')
        products = DBadd.get_products('R07618_c0')
        self.assertIn('C00579_c0', reactants) #different than below because Kbase kegg translation fie 
        self.assertIn('C00003_c0', reactants)
        self.assertIn('C15972_c0', products)
        self.assertIn('C00004_c0', products)
        self.assertIn('C00080_c0', products)

        allcpds = DBadd.get_all_compounds()
        self.assertEqual(len(allcpds), len(set(allcpds)))
        os.remove(PATH+'/kbasetestadd.db')

    def test_KEGG_with_only_KEGG_orgsinchi(self):
        """test KEGG database construction"""
        K = bkeggdb.CompileKEGGIntoDB(PATH+'/kbasetestaddinchi.db', 'bacteria', True, 4, 2, 2, 'bio', True)
        DBaddinchi = Q.Connector(PATH+'/kbasetestaddinchi.db')
        newrxns = set()
        with open('newrxns.txt') as fin:
            for line in fin:
                line = line.strip()
                newrxns.add(line)

        inter = set(newrxns) - set(kbaserxns1inchi)
        kbaserxnsadd = DBaddinchi.get_all_reactions()
        self.assertEqual(len(kbaserxns1inchi)+len(inter), len(kbaserxnsadd))

        reactants = DBaddinchi.get_reactants('R02187_c0')
        products = DBaddinchi.get_products('R02187_c0')
        self.assertIn('InChI=1S/H5O10P3/c1-11(2,3)9-13(7,8)10-12(4,5)6/h(H,7,8)(H2,1,2,3)(H2,4,5,6)_c0', reactants)
        self.assertIn('InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1_c0', reactants)
        self.assertIn('InChI=1S/H5O10P3/c1-11(2,3)9-13(7,8)10-12(4,5)6/h(H,7,8)(H2,1,2,3)(H2,4,5,6)_c0', products)
        self.assertIn('InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6-/m1/s1_c0', products)

        reactants = DBaddinchi.get_reactants('R07618_c0')
        products = DBaddinchi.get_products('R07618_c0')
        self.assertIn('C15973_c0', reactants)
        self.assertIn('InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p+1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1_c0', reactants)
        self.assertIn('C15972_c0', products)
        self.assertIn('InChI=1S/C21H29N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,25)/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1_c0', products)
        self.assertIn('InChI=1S/p+1_c0', products)

        allcpds = DBaddinchi.get_all_compounds()
        self.assertEqual(len(allcpds), len(set(allcpds)))
        os.remove(PATH+'/kbasetestaddinchi.db')

if __name__ == '__main__':
    unittest.main()
