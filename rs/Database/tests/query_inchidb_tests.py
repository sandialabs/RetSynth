from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that query the database'

import sys
import glob
import os
import re
import unittest
import cobra
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import query as Q

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)

sbml_files = glob.glob(os.path.join(PATH+'/data', '*'))
all_mets = []
all_reactions = []
all_rev = {}

'''Import model information directly from sbml files'''
print ('STATUS OF TESTS: LOADING TEST MODELS FROM SBML FILES')
for file_name in sbml_files:
    model = cobra.io.read_sbml_model(file_name)
    for r in model.reactions:
        i = r.id
        revers = r.reversibility
        all_reactions.append(i)
        if i not in all_rev:
            all_rev[i] = []
            all_rev[i].append(revers)
        elif i in all_rev and revers not in all_rev[i]:
            all_rev[i].append(revers)
    for m in model.metabolites:
        i = m.id
        all_mets.append(i)
all_mets = list(set(all_mets))
all_reactions = list(set(all_reactions))
for k, v in all_rev.iteritems():
    if len(v) > 2:
        print (' ERROR IN TEST CODE')
        sys.exit()
print ('STATUS OF TESTS: FINISHED LOADING TEST MODELS FROM SBML FILES')
init_db.Createdb(PATH+'/testinchi.db', True)
bkdb.BuildKbase(PATH+'/data', PPATH+'/KbasetoKEGGCPD.txt', PPATH+'/KbasetoKEGGRXN.txt', True, PATH+'/testinchi.db', 'bio')
DB = Q.Connector(PATH+'/testinchi.db')
os.remove(PATH+'/testinchi.db')

class Generate_databaseTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_get_all_models(self):
        print ("Testing the get_all_models() function")
        hits = DB.get_all_models()
        self.assertEqual(len(hits), 3)

    def test_get_compound_name(self):
        print ("Testing the get_compound_name() function")
        hit = DB.get_compound_name('InChI=1S/C6H10O6/c7-2-5(10)3(8)1-4(9)6(11)12/h3,5,7-8,10H,1-2H2,(H,11,12)/t3-,5+/m0/s1_c0')
        self.assertEqual('2_keto_3_deoxygluconate_c0', hit)

    def test_get_reactions(self):
        print ("Testing the get_reactions() function")
        hits = DB.get_reactions('InChI=1S/C6H10O6/c7-2-5(10)3(8)1-4(9)6(11)12/h3,5,7-8,10H,1-2H2,(H,11,12)/t3-,5+/m0/s1_c0', 1)
        self.assertIn('rxn13_c0', hits)
        hits = DB.get_reactions('InChI=1S/C6H10O6/c7-2-5(10)3(8)1-4(9)6(11)12/h3,5,7-8,10H,1-2H2,(H,11,12)/t3-,5+/m0/s1_c0', 0)
        self.assertIn('rxn7_c0', hits)

    def test_get_reaction_name(self):
        print ("Testing the get_reaction_name() function")
        hit = DB.get_reaction_name("rxn7_c0")
        self.assertEqual('rxn7_c0', hit)

    def test_get_reaction_species(self):
        print ("Testing the get_reaction_species() function")
        hit = DB.get_reaction_species('rxn7_c0')
        self.assertIn('t2', hit)

    def test_get_reactants(self):
        print ("Testing the get_reactants() function")
        hits = DB.get_reactants('rxn6_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('cpdE_c0', hits)

        hits = DB.get_reactants('rxn3_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('cpdB_c0', hits)

    def test_get_products(self):
        print ("Testing the get_products() function")
        hits = DB.get_products('rxn6_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', hits)

        hits = DB.get_products('rxn3_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('InChI=1S/2ClH.2H2N.Pt/h2*1H;2*1H2;/q;;2*-1;+4/p-2_c0', hits)

    def test_get_compounds_in_model(self):
        print ("Testing the get_compounds_in_model() function")
        hits = DB.get_compounds_in_model('t2')
        model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
        self.assertEqual(len(hits), len(model.metabolites))
        self.assertIn("cpdB_c0", hits)
        self.assertIn("InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0", hits)

    def test_get_reactions_in_model(self):
        print ("Testing the get_reactions_in_model() function")
        hits = DB.get_reactions_in_model('t2')
        model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
        self.assertEqual(len(hits), len(model.reactions))
        self.assertIn("rxn7_c0", hits)

    def test_is_reversible(self):
        print ("Testing the is_reversible() function")
        hits = DB.is_reversible('t1', "rxn3_c0")
        self.assertEqual(hits, 'true')

        hits = DB.is_reversible('t1', "rxn2_c0")
        self.assertEqual(hits, 'false')

    def test_get_all_reactions(self):
        print ("Testing the get_all_reactions() function")
        hits = DB.get_all_reactions()
        self.assertEqual(len(all_reactions), len(hits))

    def test_get_stoichiometry(self):
        print ("Testing the get_stoichiometry() function")
        hit = DB.get_stoichiometry('rxn11_c0', 'cpdI_c0', 0)
        self.assertEqual(hit[0], 1)
        hit = DB.get_stoichiometry('rxn11_c0', 'cpdH_c0', 1)
        self.assertEqual(hit[0], 1)

        hit = DB.get_stoichiometry('rxn1_c0', 'InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', 0)
        self.assertEqual(hit[0], 1)
        hit = DB.get_stoichiometry('rxn1_c0', 'cpdT_c0', 1)
        self.assertEqual(hit[0], 1)

    def test_is_reversible_all(self):
        print ("Testing the is_reversible_all() function")
        if r in all_rev:
            if len(all_rev[r]) == 2:
                hit = DB.is_reversible_all(r)
                self.assertEqual(hit, 'true')
            else:
                hit = DB.is_reversible_all(r)
                self.assertEqual(hit, all_rev[r])
        self.assertEqual(DB.is_reversible_all('rxn2_c0'), 'true')

    def test_get_compound_compartment(self):
        print ("Testing the get_compound_compartment() function")
        hit = DB.get_compound_compartment('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
        self.assertEqual(hit, 'c0')
        hit = DB.get_compound_compartment('cpdI_e0')
        self.assertEqual(hit, 'e0')

    def test_get_genes(self):
        print ("Testing the get_genes() function")
        hits = DB.get_genes('rxn11_c0', 't1')
        self.assertEqual('(g11)', hits)

        hits = DB.get_genes('rxn12_c0', 't1')
        self.assertEqual('(g2)', hits)

        hits = DB.get_genes('EX_I_e0', 't1')
        self.assertEqual(hits, 'None')

    def test_get_proteins(self):
        print ("Testing the get_proteins() function")
        hits = DB.get_proteins('rxn11_c0', 't1')
        self.assertEqual('(p11)', hits)

        hits = DB.get_proteins('rxn12_c0', 't1')
        self.assertEqual('None', hits)

        hits = DB.get_proteins('EX_I_e0', 't1')
        self.assertEqual(hits, 'None')

    def test_get_organism_name(self):
        print ("Testing the get_organism_name() function")
        hit = DB.get_organism_name('t3')
        self.assertEqual(hit, 'test3.xml')

    def test_get_uniq_metabolic_clusters(self):
        print ("Testing the get_uniq_metabolic_clusters_function()")
        hit = DB.get_uniq_metabolic_clusters()
        self.assertEqual(len(hit), 2)

    def test_get_models_from_cluster(self):
        print ("Testing the get_models_from_cluster() function")
        Q = DB.cnx.execute('SELECT DISTINCT cluster_num  FROM cluster')
        hits = list(set(Q.fetchall()))
        self.assertEqual(len(hits), 2)

    def test_get_reactions_based_on_type(self):
        hits = DB.get_reactions_based_on_type('bio')
        self.assertEqual(len(hits), 12)

    def test_get_compartment(self):
        hits = DB.get_compartment('cytosol')
        self.assertIn('c0', hits)
if __name__ == '__main__':
    unittest.main()