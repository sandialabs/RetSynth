from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests code that compares FBA results with \
                   added reactions to reaction KO FBA results'

import os
import unittest
from ShortestPath import extractinfo as ei
from FBA import build_model as bm
from FBA import optimize_target as ot
from FBA import compareKO_results as crko
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
os.remove(PATH+'/test.db')

class RuncompareKO_resultsTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.ex_info = ei.Extract_Information(osp, inmets, DB)
        self.fba = bm.BuildModel('t1', inmets, inrxns, DB)
        self.org_fbasolution = self.fba.model.optimize()
        self.optimized_fba = ot.OptimizeTarget('cpdT_c0', 't1',
                                               self.fba.model, self.ex_info.temp_rxns,
                                               self.ex_info.temp_exmets, self.fba.compounds_dict,
                                               inmets, inrxns, DB, KO=True, remove=True)
    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
        del self.ex_info
        del self.fba
        del self.optimized_fba

    def test_compareKO_results(self):
        """Test that reactions with significant flux differences are being identified"""
        print ("Test that reactions with significant flux differences are being identified")
        comparisonKOresults = crko.CompareKO('cpd08373_c0', self.optimized_fba.fbasol,
                                             self.optimized_fba.KOsolutions, self.ex_info.temp_rxns,
                                             DB)
        '''Test pathway with most flux to target compound'''
        print ("...Testing pathway with most flux to target compound")
        for r, value in comparisonKOresults.fluxchange.iteritems():
            if r in osp:
                self.assertEqual(comparisonKOresults.maxflux[r], 0)
            for rk, fluxvalue in value.iteritems():
                values = fluxvalue.split('\t')
                values[1] = float(values[1])
                values[2] = float(values[2])
                if (self.optimized_fba.fbasol.x_dict[rk] == 0 and
                        self.optimized_fba.KOsolutions[r].x_dict[rk] != 0):
                    self.assertEqual(values[1], 0)
                    self.assertTrue(values[2] > 1 or values[2] < -1)
                elif (self.optimized_fba.fbasol.x_dict[rk] != 0 and
                      self.optimized_fba.KOsolutions[r].x_dict[rk] != 0):
                    upperfold = values[1]*2.5
                    changefold = upperfold-values[1]
                    changefold = abs(changefold)
                    if values[1] < -1:
                        self.assertTrue(values[2] >= values[1]+changefold or values[2] <= upperfold)
                    elif values[1] > 1:
                        self.assertTrue(values[2] <= values[1]-changefold or values[2] >= upperfold)
if __name__ == '__main__':
    unittest.main()
