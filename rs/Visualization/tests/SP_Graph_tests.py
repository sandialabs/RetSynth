from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Test for figure generation'

import glob
import re
import os
import unittest
import shutil
from ShortestPath import extractinfo as ei
from Database import query as Q
from Visualization import SP_Graph_dot as spgd
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/Visualization/tests', '', PATH)
CPATH = re.sub('/tests', '', PATH)
if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', PPATH+'/Database/KbasetoKEGGCPD.txt',
                PPATH+'/Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
allrxns = DB.get_all_reactions()
allcpds = DB.get_all_compounds()
incpds = DB.get_compounds_in_model('t1')
inrxns = DB.get_reactions_in_model('t1')
osp = [['rxn1_c0_F', 'rxn5_c0', 'rxn6_c0', 'rxn7_c0']]
G = spgd.GraphDot(DB, PATH, incpds, inrxns)
os.remove(PATH+'/test.db')
class RunVisualizationTests(unittest.TestCase):
    def setUp(self):
        '''Initialize before every test'''
        print ("Initializing tests")
        self.ex_info = ei.Extract_Information(osp, incpds, DB)

    def tearDown(self):
        '''Clean up after each test'''
        print ("Clearing out test suite")
        del self.ex_info

    def test_graph(self):
        '''Tests assembly of nodes in a figure for a pathway'''
        print ("Tests assembly of nodes in a figure for a pathway")
        G.sc_graph('cpdT_c0', 't1', self.ex_info.temp_rxns, False)
        self.assertEqual(len(G.store_rxns)+len(G.nodes), 9)
        if os.path.isfile(PPATH+'/SC_graph_cpdT_c0_test1.xml.jpeg') is True:
            os.remove(PPATH+'/SC_graph_cpdT_c0_test1.xml.jpeg')
        if os.path.isfile(PATH+'/SC_graph_cpdT_c0_test1.xml.jpeg') is True:
            os.remove(PATH+'/SC_graph_cpdT_c0_test1.xml.jpeg')

        print ("Tests assembly of nodes in a figure for a pathway")
        G.sc_graph('cpdT_c0', 't1', self.ex_info.temp_rxns, True)
        self.assertEqual(len(G.store_rxns)+len(G.nodes), 9)
        if os.path.isfile(PPATH+'/solution_figures/SC_graph_cpdT_c0_test1.xml.jpeg') is True:
            os.remove(PPATH+'/solution_figures/SC_graph_cpdT_c0_test1.xml.jpeg')
        if os.path.isfile(PATH+'/solution_figures/SC_graph_cpdT_c0_test1.xml.jpeg') is True:
            os.remove(PATH+'/solution_figures/SC_graph_cpdT_c0_test1.xml.jpeg')
        shutil.rmtree(PATH+'/solution_figures')

        for filename in glob.glob(CPATH+"/compound*"):
            os.remove(filename)
if __name__ == '__main__':
    unittest.main()
