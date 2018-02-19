from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests constraints.py (which builds the A\
                   matrix and reaction bounds) with a test database'

import unittest
import re
import os
from ShortestPath import constraints as co
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import build_SPRESI_db as bspresidb

PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')


'''CONNECT TEST DATABASE'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', '../../Database/KbasetoKEGGCPD.txt',
                '../../Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()

init_db.Createdb(PATH+'/testRDF.db', True)
bkdb.BuildKbase(PATH+'/data', '../../Database/KbasetoKEGGCPD.txt',
                '../../Database/KbasetoKEGGRXN.txt', True,
                PATH+'/testRDF.db', 'bio')
allmets_chem = DB.get_reactions_based_on_type('chem')
bspresidb.RDF_Reader(PATH+'/data4/', PATH+'/testRDF.db', 'chem', 'c0', 1, temp_option=False, pressure_option=False,
                          yield_option=False, time_option=False, catalyst_option=False, solvent_option=False)
DBinchi = Q.Connector(PATH+'/testRDF.db')
allrxnsinchi = DBinchi.get_all_reactions()
allmetsinchi = DBinchi.get_all_compounds()
allrxns_bio = DBinchi.get_reactions_based_on_type('bio')
allrxns_chem = DBinchi.get_reactions_based_on_type('chem')
os.remove(PATH+'/test.db')
os.remove(PATH+'/testRDF.db')
print ('Not removing database')

class ConstraintsTests(unittest.TestCase):
    def setUp(self):
        '''Initialize before every test'''
        print ("Initializing tests")

    def tearDown(self):
        '''Clean up after each test'''
        print ("Clearing out test suite")

    def test_Amatrix_glpk(self):
        '''Tests A matrix construction with glpk package'''
        print ("Test A matrix construction with python glpk package")
        try:
            import glpk
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'GLPK')

            '''Test conversion of reversible reactions to single reactions'''
            print ("...Testing conversion of reversible reactions to single reactions")
            self.assertIn('rxn1_c0_F', C.allrxnsrev)
            self.assertTrue('rxn1_c0' not in C.allrxnsrev)
            self.assertEqual(len(C.lp.rows), len(C.allcpds_new))
            print (C.allrxnsrev)
            for c in C.lp.cols:
                print (c.name)
            # self.assertEqual(len(C.lp.cols), len(C.allrxnsrev))

            '''Test entries in A matrix'''
            print ("...Testing entries in A matrix")
            index_rxn1 = C.allrxnsrev.index('rxn1_c0_F')
            index_cpdT_prod = C.allcpds.index('cpdT_c0')
            index_cpdA_react = C.allcpds.index('cpdA_c0')
            index_cpdBy1_prod = C.allcpds.index('cpdBy1_c0')
            index_cpdC_react = C.allcpds.index('cpdC_c0')
            index_cpdB_prod = C.allcpds.index('cpdB_c0')
            index_rxn3 = C.allrxnsrev.index('rxn3_c0')

            self.assertEqual(C.A[index_cpdA_react][index_rxn1], -1)
            self.assertEqual(C.A[index_cpdT_prod][index_rxn1], 1)
            self.assertEqual(C.A[index_cpdBy1_prod][index_rxn1], 1)

            self.assertEqual(C.A[index_cpdB_prod][index_rxn3], 1)
            self.assertEqual(C.A[index_cpdC_react][index_rxn3], -1)
        except ImportError:
            print ('pyglpk package not instaled, cannot run tests for this python package')

    def test_Amatrix_glpk_bio(self):
        '''Tests A matrix construction with glpk package'''
        print ("Test A matrix construction with python glpk package")
        try:
            import glpk
            C = co.ConstructInitialLP(allrxns_bio, allmetsinchi, DBinchi, [], True, False, 'GLPK')

            '''Test conversion of reversible reactions to single reactions'''
            print ("...Testing conversion of reversible reactions to single reactions")
            self.assertIn('rxn1_c0_F', C.allrxnsrev)
            self.assertTrue('rxn1_c0' not in C.allrxnsrev)
            self.assertEqual(len(C.lp.rows), len(C.allcpds_new))
            self.assertEqual(len(C.lp.cols), len(C.allrxnsrev))

            '''Test entries in A matrix'''
            print ("...Testing entries in A matrix")
            index_rxn1 = C.allrxnsrev.index('rxn1_c0_F')
            index_cpdT_prod = C.allcpds.index('cpdT_c0')
            index_cpdA_react = C.allcpds.index('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
            index_cpdBy1_prod = C.allcpds.index('cpdBy1_c0')
            index_cpdC_react = C.allcpds.index('InChI=1S/2ClH.2H2N.Pt/h2*1H;2*1H2;/q;;2*-1;+4/p-2_c0')
            index_cpdB_prod = C.allcpds.index('cpdB_c0')
            index_rxn3 = C.allrxnsrev.index('rxn3_c0')

            self.assertEqual(C.A[index_cpdA_react][index_rxn1], -1)
            self.assertEqual(C.A[index_cpdT_prod][index_rxn1], 1)
            self.assertEqual(C.A[index_cpdBy1_prod][index_rxn1], 1)

            self.assertEqual(C.A[index_cpdB_prod][index_rxn3], 1)
            self.assertEqual(C.A[index_cpdC_react][index_rxn3], -1)
        except ImportError:
            print ('pyglpk package not instaled, cannot run tests for this python package')

    def test_Amatrix_glpk_chem(self):
        '''Tests A matrix construction with glpk package'''
        print ("Test A matrix construction with python glpk package")
        try:
            import glpk
            C = co.ConstructInitialLP(allrxns_chem, allmetsinchi, DBinchi, [], True, False, 'GLPK')
            self.assertEqual(len(allrxns_chem), 9)
            '''Test conversion of reversible reactions to single reactions'''
            print ("...Testing conversion of reversible reactions to single reactions")
            self.assertIn('rxn4073836_s', C.allrxnsrev)
            self.assertEqual(len(C.lp.rows), len(C.allcpds_new))
            self.assertEqual(len(C.lp.cols), len(C.allrxnsrev))

            '''Test entries in A matrix'''
            print ("...Testing entries in A matrix")
            index_rxn = C.allrxnsrev.index('rxn4073836_s')
            index_react = C.allcpds.index('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0')
            index_prod = C.allcpds.index('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0')

            self.assertEqual(C.A[index_react][index_rxn], -1)
            self.assertEqual(C.A[index_prod][index_rxn], 1)

        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')


    def test_Amatrix_pulp_bio(self):
        '''Tests A matrix construction with glpk package'''
        print ("Test A matrix construction with python glpk package")
        try:
            import pulp
            C = co.ConstructInitialLP(allrxns_bio, allmetsinchi, DBinchi, [], True, False, 'PULP')

            '''Test conversion of reversible reactions to single reactions'''
            print ("...Testing conversion of reversible reactions to single reactions")
            self.assertIn('rxn1_c0_F', C.allrxnsrev)
            self.assertTrue('rxn1_c0' not in C.allrxnsrev)
            self.assertEqual(len(C.lp.constraints), len(C.allcpds_new))
            self.assertEqual(len(C.variables), len(C.allrxnsrev))

            '''Test entries in A matrix'''
            print ("...Testing entries in A matrix")
            index_rxn1 = C.allrxnsrev.index('rxn1_c0_F')
            index_cpdT_prod = C.allcpds.index('cpdT_c0')
            index_cpdA_react = C.allcpds.index('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
            index_cpdBy1_prod = C.allcpds.index('cpdBy1_c0')
            index_cpdC_react = C.allcpds.index('InChI=1S/2ClH.2H2N.Pt/h2*1H;2*1H2;/q;;2*-1;+4/p-2_c0')
            index_cpdB_prod = C.allcpds.index('cpdB_c0')
            index_rxn3 = C.allrxnsrev.index('rxn3_c0')

            self.assertEqual(C.A[index_cpdA_react][index_rxn1], -1)
            self.assertEqual(C.A[index_cpdT_prod][index_rxn1], 1)
            self.assertEqual(C.A[index_cpdBy1_prod][index_rxn1], 1)

            self.assertEqual(C.A[index_cpdB_prod][index_rxn3], 1)
            self.assertEqual(C.A[index_cpdC_react][index_rxn3], -1)
        except ImportError:
            print ('pyglpk package not installed, cannot run tests for this python package')

    def test_Amatrix_pulp_chem(self):
        '''Tests A matrix construction with glpk package'''
        print ("Test A matrix construction with python glpk package")
        try:
            import pulp
            C = co.ConstructInitialLP(allrxns_chem, allmetsinchi, DBinchi, [], True, False, 'PULP')
            self.assertEqual(len(allrxns_chem), 9)
            '''Test conversion of reversible reactions to single reactions'''
            print ("...Testing conversion of reversible reactions to single reactions")
            self.assertIn('rxn4073836_s', C.allrxnsrev)
            self.assertEqual(len(C.lp.constraints), len(C.allcpds_new))
            self.assertEqual(len(C.variables), len(C.allrxnsrev))

            '''Test entries in A matrix'''
            print ("...Testing entries in A matrix")
            index_rxn = C.allrxnsrev.index('rxn4073836_s')
            index_react = C.allcpds.index('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0')
            index_prod = C.allcpds.index('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0')

            self.assertEqual(C.A[index_react][index_rxn], -1)
            self.assertEqual(C.A[index_prod][index_rxn], 1)

        except ImportError:
            print ('pyglpk package not instaled, cannot run tests for this python package')

    def test_Amatrix_pulp(self):
        '''Tests A matrix construction with pulp package'''
        print ("Test A matrix construction with python pulp package")
        try:
            import pulp
            C = co.ConstructInitialLP(allrxns, allmets, DB, [], True, False, 'PULP')

            '''Test conversion of reversible reactions to single reactions'''
            print ("...Testing conversion of reversible reactions to single reactions")
            self.assertIn('rxn1_c0_F', C.allrxnsrev)
            self.assertTrue('rxn1_c0' not in C.allrxnsrev)
            self.assertEqual(len(C.lp.constraints), len(C.allcpds_new))
            self.assertEqual(len(C.variables), len(C.allrxnsrev))

            '''Test entries in A matrix'''
            print ("...Testing entries in A matrix")
            index_rxn1 = C.allrxnsrev.index('rxn1_c0_F')
            index_cpdT_prod = C.allcpds.index('cpdT_c0')
            index_cpdA_react = C.allcpds.index('cpdA_c0')
            index_cpdBy1_prod = C.allcpds.index('cpdBy1_c0')
            index_cpdC_react = C.allcpds.index('cpdC_c0')
            index_cpdB_prod = C.allcpds.index('cpdB_c0')
            index_rxn3 = C.allrxnsrev.index('rxn3_c0')

            self.assertEqual(C.A[index_cpdA_react][index_rxn1], -1)
            self.assertEqual(C.A[index_cpdT_prod][index_rxn1], 1)
            self.assertEqual(C.A[index_cpdBy1_prod][index_rxn1], 1)

            self.assertEqual(C.A[index_cpdB_prod][index_rxn3], 1)
            self.assertEqual(C.A[index_cpdC_react][index_rxn3], -1)
        except ImportError:
            print ('pulp package not instaled, cannot run tests for this python package')

if __name__ == '__main__':
    unittest.main()
