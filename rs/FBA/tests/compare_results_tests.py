from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests code that compares wild type (no added reactions) \
                   FBA results to FBA results with added reactions'

import os
import unittest
import sys
import re
from ShortestPath import extractinfo as ei
from FBA import build_model as bm
from FBA import optimize_target as ot
from FBA import compare_results as cr
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

'''Load database'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', '../../Database/KbasetoKEGGCPD.txt',
                '../../Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')

allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
inmets = DB.get_compounds_in_model('t1')
inrxns = DB.get_reactions_in_model('t1')
osp = [['rxn1_c0_F', 'rxn2_c0', 'rxn3_c0_F']]
ex_info = ei.Extract_Information(osp, inmets, DB)
fba = bm.BuildModel('t1', inmets, inrxns, DB)
org_fbasolution = fba.model.optimize()
optimized_fba = ot.OptimizeTarget('cpdT_c0', 't1', fba.model,
                                  ex_info.temp_rxns, ex_info.temp_exmets,
                                  fba.compounds_dict, inmets, inrxns, DB,
                                  KO=False, remove=True)
os.remove(PATH+'/test.db')
class Runcompare_resultsTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
    def test_compare_results(self):
        """Test that reactions with significant flux differences are being identified"""
        print ("Testing that reactions with significant flux differences are being identified")
        comparisonresults = cr.Compare('cpdT_c0', org_fbasolution, optimized_fba.fbasol,
                                       ex_info.temp_rxns, DB)
        '''Test pathway with most flux to target compound'''
        print ("...Testing pathway with most flux to target compound")
        self.assertEqual(comparisonresults.maxpath, 1)

        '''Ensure that flux is actually flowing through pathway'''
        print ("...Testing that flux is actually flowing through pathway")
        self.assertTrue(comparisonresults.maxflux > 0)
        for r, fluxvalue in comparisonresults.fluxchange.iteritems():
            values = fluxvalue.split('\t')
            match = re.search('$\w+', values[0])
            if match is not None:
                values[0] = float(values[0])
                values[1] = float(values[1])
                if org_fbasolution.x_dict[r] == 0 and optimized_fba.model.solution.x_dict[r] != 0:
                    self.assertEqual(values[0], 0)
                    self.assertTrue(values[1] > 1 or values[1] < -1)
                elif (self.org_fbasolution.x_dict[r] != 0 and
                      optimized_fba.model.solution.x_dict[r] != 0):
                    upperfold = values[0]*2.5
                    changefold = upperfold-values[0]
                    changefold = abs(changefold)
                    if values[0] < 0:
                        self.assertTrue(values[1] >= values[0]+changefold or
                                        values[1] <= upperfold)
                    elif values[0] > 0:
                        self.assertTrue(values[1] <= values[0]-changefold or
                                        values[1] >= upperfold)
                else:
                    print ('Test should be fail')
                    sys.exit()
            else:
                if r not in org_fbasolution.x_dict:
                    self.assertTrue(values[0] != 0)
if __name__ == '__main__':
    unittest.main()
