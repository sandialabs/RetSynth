from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests optimizing of a FBA model from the \
                   Database and addition of external reactions to model'

import os
import re
import unittest
import cobra
from FBA import build_model as bm
from FBA import optimize_target as ot
from Database import query as Q
from ShortestPath import extractinfo as ei
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/FBA/tests', '', PATH)

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')
'''Load database'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data1', PPATH+'/Database/data/KbasetoKEGGCPD.txt',
                PPATH+'/Database/data/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')

allrxns = DB.get_all_reactions()
allmets = DB.get_all_compounds()
inmets = DB.get_compounds_in_model('t1')
inrxns = DB.get_reactions_in_model('t1')

model = cobra.io.read_sbml_model(PATH+'/data1/test1.xml')

class RunOptimizeTargetTests(unittest.TestCase):
    def setUp(self):
        """Initalizes tests"""
        print ("Initializing tests")

    def tearDown(self):
        """clears out tests"""
        print ("Clearing out test suite")

    def test_optimize_target(self):
        '''Testing the optimization of target compound using pulp'''
        print ('Testing the optimization of target compound using pulp')
        optimal_pathways = [['rxn1_c0_F', 'rxn2_c0_R', 'rxn3_c0_R']]
        ex_info = ei.Extract_Information(optimal_pathways, inmets, inrxns, DB)
        fba = bm.BuildModel('t1', inmets, inrxns, DB, False)

        print ("Testing the FBA optimization of a target compound")
        self.assertEqual(len(model.metabolites), len(fba.model.metabolites))
        self.assertEqual(len(model.reactions), len(fba.model.reactions))

        optimized_fba = ot.OptimizeTarget('cpdT_c0', 't1', fba.model,
                                          ex_info.temp_rxns, ex_info.temp_exmets,
                                          fba.compounds_dict, inmets, inrxns, DB,
                                          False, KO=True, remove=False)

        print ("...Test correct number of reactions were added'")
        self.assertEqual(len(optimized_fba.model.reactions),
                         (len(model.reactions)+len(ex_info.temp_rxns[1])+
                          len(optimized_fba.sink_rxns)))

        print ("...Testing that correct number of metabolites were added")
        self.assertEqual(len(optimized_fba.model.metabolites),
                         (len(model.metabolites)+len(ex_info.temp_exmets[1])))

        print ("...Testing knockouts")
        self.assertTrue(len(optimized_fba.KOsolutions) > 0)

        print ("...Testing that essential reactions were identified")
        self.assertTrue(len(optimized_fba.essentialrxns) > 0)
        self.assertIn('rxn1_c0', optimized_fba.essentialrxns)
if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/test.db')
