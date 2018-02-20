from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes add reactions and compounds to existing database'

import os
import re
import sqlite3
import unittest
from Database import query as Q
from Database import initialize_database as init_db
from Database.build_metacyc_db import MetaCyc
from Database import build_metacyc_db as bmcdb
from Database import build_kbase_db as bkdb

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)

if os.path.isfile(PATH+'/kbasetest.db') is True:
    os.remove(PATH+'/kbasetest.db')
init_db.Createdb(PATH+'/kbasetest.db', False)
bkdb.BuildKbase(PATH+'/datam', PPATH+'/KbasetoKEGGCPD.txt', PPATH+'/KbasetoKEGGRXN.txt', False, PATH+'/kbasetest.db', 'bio')
DB = Q.Connector(PATH+'/kbasetest.db')
file_name = open(PATH+'/metacyc_data/MetaCyc.aliases')
line = file_name.readline()
BIOCYC_translator = {}
BIOCYC_translator['rxn'] = {}
BIOCYC_translator['compound'] = {}
if line.startswith('#'):
    pass
for count, line in enumerate(file_name):
    larray = line.strip('\n').split('\t')
    larray[0] = re.sub('\.\w+$', '', larray[0])
    if larray[1] != '':
        if larray[0] not in BIOCYC_translator['rxn'].keys():
            BIOCYC_translator['rxn'][larray[0]] = []
            BIOCYC_translator['rxn'][larray[0]].append(larray[1])
        else:
            BIOCYC_translator['rxn'][larray[0]].append(larray[1])
    elif larray[2] != '':
        if larray[0] not in BIOCYC_translator['compound']:
            BIOCYC_translator['compound'][larray[0]] = []
            BIOCYC_translator['compound'][larray[0]].append(larray[2])
        else:
            BIOCYC_translator['compound'][larray[0]].append(larray[2])

class Translate_metacycTests(unittest.TestCase):
    def test_metacyc(self):
        """Tests that metacyc can be added to a kbase database correctly"""
        print ("Testing that metacyc can be added to a kbase database correctly")
        MC = MetaCyc(DB, False, sqlite3.connect(PATH+'/kbasetest.db'))
        MC.read_metacyc_file(BIOCYC_translator, PATH+'/metacyc_data/metabolic-reactions.xml')
        self.assertIn(('R04139_c0', 'META', '(G-44401)'), MC.genelist)
        self.assertIn(('R04139_c0', 'META', 'true'), MC.modelreactions)
        self.assertIn(('R00850_c0', 'META', 'true'), MC.modelreactions)
        self.assertIn(('R00850_c0_v1', 'META', 'true'), MC.modelreactions)
        self.assertIn(('R00850_c0_v2', 'META', 'true'), MC.modelreactions)
        self.assertIn(('R00850_c0_v3', 'META', 'true'), MC.modelreactions)
        self.assertEqual(MC.all_rxnreversibility['R04139_c0'], 'true')
        self.assertEqual(MC.all_rxnreversibility['R00850_c0'], 'true')
        self.assertEqual(MC.all_rxnreversibility['R00850_c0_v1'], 'true')
        self.assertEqual(MC.all_rxnreversibility['R00850_c0_v2'], 'true')
        self.assertEqual(MC.all_rxnreversibility['R00850_c0_v3'], 'true')
        self.assertEqual(5, len(MC.all_reaction_compound['R04139_c0']))

        self.assertIn(('R04139_c0', 'C00003_c0', False, '1'),
                      MC.all_reaction_compound['R04139_c0'])
        self.assertIn(('R04139_c0', 'C00004_c0', True, '1'),
                      MC.all_reaction_compound['R04139_c0'])
        self.assertIn(('R04139_c0', 'C00003_c0', False, '1'),
                      MC.all_reaction_compound['R04139_c0'])
        self.assertIn(('TRANS__45__RXN0__45__443_t0', 'C00638_p0', False, '1'),
                      MC.all_reaction_compound['TRANS__45__RXN0__45__443_t0'])
        self.assertIn(('TRANS__45__RXN0__45__443_t0', 'C00638_c0', True, '1'),
                      MC.all_reaction_compound['TRANS__45__RXN0__45__443_t0'])

        self.assertIn(('R06683_c0', 'C12424_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00005_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00080_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00007_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C12425_c0', True, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00006_c0', True, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C12425_c0', True, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'META', '(G-17686)'),
                      MC.genelist)
        self.assertIn(('R01173_c0', 'META', '(CA_P0035) or (CA_P0162)'), MC.genelist)
        self.assertIn(('R01173_c0_v1', 'META', '(CA_P0035) or (CA_P0162)'), MC.genelist)
        self.assertIn(('R06683_c0', 'META', '1.14.13.180'),
                      MC.proteinlist)
        self.assertIn(('R00850_c0', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn(('R00850_c0_v1', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn(('R00850_c0_v2', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn(('R00850_c0_v3', 'META', '2.7.1.142'), MC.proteinlist)
        os.remove(PATH+'/kbasetest.db')

    def test_BIOCYC_translator(self):
        """Test that translator dictionary is built correctly"""
        print ("Testing that translator dictionary is built correctly")
        T = bmcdb.Translate(PATH+'/kbasetest.db', DB, PATH+'/metacyc_data/metabolic-reactions.xml',
                         PATH+'/metacyc_data/MetaCyc.aliases', False, 'bio', add=False)
        self.assertEqual(len(T.BIOCYC_translator['rxn']), len(BIOCYC_translator['rxn']), 'bio')
        self.assertEqual(len(T.BIOCYC_translator['compound']), len(BIOCYC_translator['compound']))
        self.assertEqual(len(BIOCYC_translator['compound'])+len(BIOCYC_translator['rxn']), 26233)
        self.assertEqual(T.BIOCYC_translator['compound']['1-RADYL-2-ACYL-SN-GLYCERO-3-PHOSPHOLIPID'],
                         ['cpd21786'])
        self.assertIn('rxn16863', T.BIOCYC_translator['rxn']['1.1.1.190-RXN'])
        self.assertIn('rxn01938', T.BIOCYC_translator['rxn']['1.1.1.190-RXN'])
        self.assertEqual(len(T.BIOCYC_translator['rxn']['1.1.1.190-RXN']), 2)

if __name__ == '__main__':
    unittest.main()
