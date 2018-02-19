from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests read_targets code'

import os
import unittest
from Parser import read_targets as rt
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb

PATH = os.path.dirname(os.path.abspath(__file__))

'''Load database'''
if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')
if os.path.isfile(PATH+'/testinchi.db') is True:
    os.remove(PATH+'/testinchi.db')

init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', '../../Database/KbasetoKEGGCPD.txt',
                '../../Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')

init_db.Createdb(PATH+'/testinchi.db', True)
bkdb.BuildKbase(PATH+'/data', '../../Database/KbasetoKEGGCPD.txt',
                '../../Database/KbasetoKEGGRXN.txt', True,
                PATH+'/testinchi.db', 'bio')
DBinchi = Q.Connector(PATH+'/testinchi.db')


answer1 = [['cpdT_c0', '', 't1', '']]
answer2 = [['cpdA_c0', '71089297', 't1', '']]
answer3 = [['cpdA_c0', '71089297', 't1', 'test1']]
answer4 = [['cpdT_c0', 'cpdT_c0', 't1', 'test1']]
answer5 = [['cpdT_c0', 'cpdT_c0', 't1', ''],
           ['cpdT_c0', 'cpdT_c0', 't3', '']]

answer1inchi = [['InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', '', 't1', '']]
answer2inchi = [['InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', '71089297', 't1', '']]
answer3inchi = [['InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', '71089297', 't1', 'test1']]
answer4inchi = [['InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', 'cpdA_c0', 't1', 'test1']]
answer5inchi = [['InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', 'cpdA_c0', 't1', ''],
                ['InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', 'cpdA_c0', 't3', '']]

class PubchemTests(unittest.TestCase):
    def test_read_targets1(self):
        """Test targets inputed into rs"""
        R = rt.Readfile(PATH+'/data2/test_targets1.txt', DB)
        print ("Test targets inputed into rs")
        self.assertEqual(R.targets, answer1)

    def test_read_targets2(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data2/test_targets2.txt', DB)
        self.assertEqual(R.targets, answer2)

    def test_read_targets3(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data2/test_targets3.txt', DB)
        self.assertEqual(R.targets, answer3)

    def test_read_targets4(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data2/test_targets4.txt', DB)
        self.assertEqual(R.targets, answer4)

    def test_read_targets5(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data2/test_targets5.txt', DB)
        self.assertEqual(R.targets, answer5)

    def test_read_targets1inchi(self):
        """Test targets inputed into rs"""
        R = rt.Readfile(PATH+'/data3/test_targets1.txt', DBinchi)
        print ("Test targets inputed into rs")
        self.assertEqual(R.targets, answer1inchi)

    def test_read_targets2inchi(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data3/test_targets2.txt', DBinchi)
        self.assertEqual(R.targets, answer2inchi)

    def test_read_targets3inchi(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data3/test_targets3.txt', DBinchi)
        self.assertEqual(R.targets, answer3inchi)

    def test_read_targets4inchi(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data3/test_targets4.txt', DBinchi)
        self.assertEqual(R.targets, answer4inchi)

    def test_read_targets5inchi(self):
        """Test targets inputed into rs"""
        print ("Test targets inputed into rs")
        R = rt.Readfile(PATH+'/data3/test_targets5.txt', DBinchi)
        self.assertEqual(R.targets, answer5inchi)
        os.remove(PATH+'/testinchi.db')
        os.remove(PATH+'/test.db')
if __name__ == '__main__':
    unittest.main()
