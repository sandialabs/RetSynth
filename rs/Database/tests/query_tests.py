from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that query the database'

import sys
import glob
import os
import unittest
import cobra
from Database import query as Q
from Database import generate_database as gen_db

PATH = os.path.dirname(os.path.abspath(__file__))

'''GENERATE TEST DATABASE'''
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
gen_db.Createdb(PATH+'/test.db', PATH+'/data', False, 'bio')

class Generate_databaseTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.DB = Q.Connector(PATH+'/test.db')

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
        del self.DB

    def test_get_all_models(self):
        print ("Testing the get_all_models() function")
        hits = self.DB.get_all_models()
        self.assertEqual(len(hits), 3)

    def test_get_compound_name(self):
        print ("Testing the get_compound_name() function")
        hit = self.DB.get_compound_name('cpdF_c0')
        self.assertEqual('2_keto_3_deoxygluconate_c0', hit)

    def test_get_reactions(self):
        print ("Testing the get_reactions() function")
        hits = self.DB.get_reactions('cpdF_c0', 1)
        self.assertIn('rxn13_c0', hits)
        hits = self.DB.get_reactions('cpdF_c0', 0)
        self.assertIn('rxn7_c0', hits)

    def test_get_reaction_name(self):
        print ("Testing the get_reaction_name() function")
        hit = self.DB.get_reaction_name("rxn7_c0")
        self.assertEqual('rxn7_c0', hit)

    def test_get_reaction_species(self):
        print ("Testing the get_reaction_species() function")
        hit = self.DB.get_reaction_species('rxn7_c0')
        self.assertIn('t2', hit)

    def test_get_reactants(self):
        print ("Testing the get_reactants() function")
        hits = self.DB.get_reactants('rxn6_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('cpdE_c0', hits)

    def test_get_products(self):
        print ("Testing the get_products() function")
        hits = self.DB.get_products('rxn6_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('cpdD_c0', hits)

    def test_get_compounds_in_model(self):
        print ("Testing the get_compounds_in_model() function")
        hits = self.DB.get_compounds_in_model('t2')
        model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
        self.assertEqual(len(hits), len(model.metabolites))
        self.assertIn("cpdB_c0", hits)

    def test_get_reactions_in_model(self):
        print ("Testing the get_reactions_in_model() function")
        hits = self.DB.get_reactions_in_model('t2')
        model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
        self.assertEqual(len(hits), len(model.reactions))
        self.assertIn("rxn7_c0", hits)

    def test_is_reversible(self):
        print ("Testing the is_reversible() function")
        hits = self.DB.is_reversible('t1', "rxn3_c0")
        self.assertEqual(hits, 'true')

        hits = self.DB.is_reversible('t1', "rxn2_c0")
        self.assertEqual(hits, 'false')

    def test_get_all_compounds(self):
        print ("Testing the get_all_compounds() function")
        hits = self.DB.get_all_compounds()
        self.assertEqual(len(hits), len(all_mets))

    def test_get_all_reactions(self):
        print ("Testing the get_all_reactions() function")
        hits = self.DB.get_all_reactions()
        self.assertEqual(len(all_reactions), len(hits))

    def test_get_stoichiometry(self):
        print ("Testing the get_stoichiometry() function")
        hit = self.DB.get_stoichiometry('rxn11_c0', 'cpdI_c0', 0)
        self.assertEqual(hit[0], 1)
        hit = self.DB.get_stoichiometry('rxn11_c0', 'cpdH_c0', 1)
        self.assertEqual(hit[0], 1)

    def test_is_reversible_all(self):
        print ("Testing the is_reversible_all() function")
        if r in all_rev:
            if len(all_rev[r]) == 2:
                hit = self.DB.is_reversible_all(r)
                self.assertEqual(hit, 'true')
            else:
                hit = self.DB.is_reversible_all(r)
                self.assertEqual(hit, all_rev[r])
        self.assertEqual(self.DB.is_reversible_all('rxn2_c0'), 'true')
        os.remove(PATH+'/test.db')

    def test_get_compound_compartment(self):
        print ("Testing the get_compound_compartment() function")
        hit = self.DB.get_compound_compartment('cpdA_c0')
        self.assertEqual(hit, 'c0')
        hit = self.DB.get_compound_compartment('cpdI_e0')
        self.assertEqual(hit, 'e0')

    def test_get_genes(self):
        print ("Testing the get_genes() function")
        hits = self.DB.get_genes('rxn11_c0', 't1')
        self.assertEqual('(g11)', hits)

        hits = self.DB.get_genes('rxn12_c0', 't1')
        self.assertEqual('(g2)', hits)

        hits = self.DB.get_genes('EX_I_e0', 't1')
        self.assertEqual(hits, 'None')

    def test_get_proteins(self):
        print ("Testing the get_proteins() function")
        hits = self.DB.get_proteins('rxn11_c0', 't1')
        self.assertEqual('(p11)', hits)

        hits = self.DB.get_proteins('rxn12_c0', 't1')
        self.assertEqual('None', hits)

        hits = self.DB.get_proteins('EX_I_e0', 't1')
        self.assertEqual(hits, 'None')

    def test_get_organism_name(self):
        print ("Testing the get_organism_name() function")
        hit = self.DB.get_organism_name('t3')
        self.assertEqual(hit, 'test3.xml')

    def test_get_uniq_metabolic_clusters(self):
        print ("Testing the get_uniq_metabolic_clusters_function()")
        hit = self.DB.get_uniq_metabolic_clusters()
        self.assertEqual(len(hit), 2)

    def test_get_models_from_cluster(self):
        print ("Testing the get_models_from_cluster() function")
        Q = self.DB.cnx.execute('SELECT DISTINCT cluster_num  FROM cluster')
        hits = list(set(Q.fetchall()))
        self.assertEqual(len(hits), 2)

    def test_get_reactions_based_on_type(self):
        hits = self.DB.get_reactions_based_on_type('bio')
        self.assertEqual(len(hits), 12)

    def test_get_compartment(self):
        hits = self.DB.get_compartment('cytosol')
        self.assertIn('c0', hits)
if __name__ == '__main__':
    unittest.main()
