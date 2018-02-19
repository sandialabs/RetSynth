from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that generate database on test set of metabolic networks'

import glob
import os
import sys
import sqlite3
import unittest
from Database import build_kbase_db as bkdb
import cobra
'''Build test database'''
PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')
sqlite_database = sqlite3.connect(PATH+'/test.db')
sqlite_database.execute('''CREATE table model (ID text,file_name text)''')
sqlite_database.execute('''CREATE table compound (ID text, name text,compartment text, kegg_id)''')
sqlite_database.execute('''CREATE table model_compound (cpd_ID text,model_ID)''')
sqlite_database.execute('''CREATE table reaction (ID text, name text, kegg_id, type text)''')
sqlite_database.execute('''CREATE table compartments (ID text, name text)''')

sqlite_database.execute('''CREATE table model_reaction
                       (reaction_ID text, model_ID text, is_rev bit(1))''')
sqlite_database.execute('''CREATE table reaction_compound
                       (reaction_ID text, cpd_ID text, is_prod bit(1),
                        stoichiometry int, filenum int)''')
sqlite_database.execute('''CREATE table reaction_reversibility
                       (reaction_ID text, is_reversible bit(1))''')
sqlite_database.execute('''CREATE table reaction_gene
                       (reaction_ID text,model_ID text,gene_ID text)''')
sqlite_database.execute('''CREATE table reaction_protein
                       (reaction_ID text, model_ID text, protein_ID text)''')
sqlite_database.execute('''CREATE table cluster (cluster_num text, ID text)''')
sqlite_database.execute('''CREATE table original_db_cpdIDs (ID text, inchi_id text)''')

sqlite_database.execute('''CREATE INDEX reactioncompound_ind1 ON
                        reaction_compound(reaction_ID, cpd_ID, is_prod)''')
sqlite_database.execute('''CREATE INDEX reactioncompound_ind2 ON
                        reaction_compound(cpd_ID, is_prod)''')
sqlite_database.execute('''CREATE INDEX modelreaction_ind1 ON
                        model_reaction(model_ID)''')
sqlite_database.execute('''CREATE INDEX modelreaction_ind2 ON
                        model_reaction(reaction_ID)''')
sqlite_database.execute('''CREATE INDEX modelcompound_ind1 ON
                        model_compound(model_ID)''')
sqlite_database.execute('''CREATE INDEX modelcompound_ind2 ON
                        model_compound(cpd_ID)''')
sqlite_database.execute('''CREATE INDEX model_ind ON model(ID)''')
sqlite_database.execute('''CREATE INDEX reaction_ind ON reaction(ID)''')
sqlite_database.execute('''CREATE INDEX reaction_reversibility_ind ON
                        reaction_reversibility(reaction_ID)''')
sqlite_database.execute('''CREATE INDEX compound_ind ON compound(ID)''')
sqlite_database.execute('''CREATE INDEX reactiongene_ind ON
                        reaction_gene(reaction_ID, model_ID)''')
sqlite_database.execute('''CREATE INDEX reactionprotein_ind ON
                        reaction_protein(reaction_ID,model_ID)''')
sqlite_database.execute('''CREATE INDEX cluster_ind ON cluster(cluster_num)''')
sqlite_database.execute('''CREATE INDEX original_db_cpdIDs_ind ON
                        original_db_cpdIDs(ID, inchi_id)''')
sbml_files = glob.glob(os.path.join(PATH+'/data', '*'))
all_mets = []
all_reactions = []
all_rev = {}

'''Import model information directly from sbml files'''
print ("""STATUS OF TESTS: Loading test models from SBML files using
      CobraPy to compare to information in the database""")

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
print ('STATUS OF TESTS: Finished loading test models from SBML files using CobraPy')

class Generate_databaseTests_PubchemTrue(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_build_kbase(self):
        ''' test construction of kbase tables'''
        bkdb.BuildKbase(PATH+'/data', '../KbasetoKEGGCPD.txt', '../KbasetoKEGGRXN.txt', False, PATH+'/test.db', 'bio')
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute('SELECT * FROM model')
        results_model = Q.fetchall()
        self.assertEqual(len(results_model), len(sbml_files))

        print ("""Testing to see if database is filling tables with correct number of
                      reactions, compounds, reversibility, and reaction compound information""")
        Q = cnx.execute('SELECT * FROM reaction')
        results_rxns = Q.fetchall()
        self.assertEqual(len(results_rxns), len(all_reactions))
        value = ('rxn1_c0', 'rxn1_c0', 'None', 'bio')
        value1 = ('rxn2_c0', 'rxn2_c0', 'None', 'bio')
        value2 = ('rxn11_c0', 'rxn11_c0', 'None', 'bio')
        self.assertIn(value, results_rxns)
        self.assertIn(value1, results_rxns)
        self.assertIn(value2, results_rxns)

        print ("...Testing clustering")
        Q = cnx.execute('SELECT * FROM cluster')
        results_cluster = Q.fetchall()
        cluster_1_num = []
        for result in results_cluster:
            if result[1] == 't3' or result[1] == 't1':
                cluster_1_num.append(result[0])
        Q = cnx.execute('SELECT DISTINCT cluster_num FROM cluster')
        results_cluster_unique = Q.fetchall()
        self.assertEqual(cluster_1_num[0], cluster_1_num[1])
        self.assertEqual(len(results_cluster_unique), 2)
        self.assertEqual(len(results_cluster), 3)

        print ("...Testing total number of metabolites inputed into database")
        Q = cnx.execute('SELECT * FROM compound')
        results_compound = Q.fetchall()
        self.assertEqual(len(results_compound), len(all_mets))
        self.assertIn(('cpdA_c0', 'cpdA_c0', 'c0', 'None'), results_compound)
        self.assertIn(('cpdB_c0', 'cpdB_c0', 'c0', 'None',), results_compound)
        self.assertIn(('cpdF_c0', '2_keto_3_deoxygluconate_c0', 'c0', 'None'), results_compound)
        self.assertIn(('cpdI_e0', 'cpdI_e0', 'e0',  'None'), results_compound)

        print ("...Testing total number of reaction reversibility records inputed into database")
        Q = cnx.execute('SELECT * FROM reaction_reversibility')
        results_rxns_revers = Q.fetchall()
        self.assertEqual(len(results_rxns_revers), len(all_reactions))
        self.assertIn(('rxn1_c0', 'true'), results_rxns_revers)
        self.assertIn(('rxn2_c0', 'true'), results_rxns_revers)
        self.assertIn(('rxn3_c0', 'true'), results_rxns_revers)
        self.assertIn(('rxn4_c0', 'false'), results_rxns_revers)

        print ("...Testing reaction compound relationship information inputed into database")
        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID='rxn2_c0'")
        results_cpd = Q.fetchall()
        results_cpd2 = [tuple(y for c, y in enumerate(x) if c != 4) for x in results_cpd]
        value = ('rxn2_c0', 'cpdA_c0', 0, 1)
        value1 = ('rxn2_c0', 'cpdB_c0', 1, 1)
        value2 = ('rxn2_c0', 'cpdBy1_c0', 1, 1)

        valuev2 = ('rxn2_c0', 'cpdA_c0', 1, 1)
        value1v2 = ('rxn2_c0', 'cpdB_c0', 0, 1)
        value2v2 = ('rxn2_c0', 'cpdBy1_c0', 0, 1)
        self.assertEqual(len(results_cpd), 3)

        try:
            self.assertIn(value, results_cpd2)
            self.assertIn(value1, results_cpd2)
            self.assertIn(value2, results_cpd2)

        except AssertionError:
            self.assertIn(valuev2, results_cpd2)
            self.assertIn(value1v2, results_cpd2)
            self.assertIn(value2v2, results_cpd2)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID='rxn3_c0'")
        results_cpd = Q.fetchall()
        results_cpd2 = [tuple(y for c, y in enumerate(x) if c != 4) for x in results_cpd]
        value = ('rxn3_c0', 'cpdB_c0', 0, 1)
        value1 = ('rxn3_c0', 'cpdC_c0', 1, 1)

        valuev2 = ('rxn3_c0', 'cpdB_c0', 1, 1)
        value1v2 = ('rxn3_c0', 'cpdC_c0', 0, 1)

        self.assertEqual(len(results_cpd), 2)
        try:
            self.assertIn(value, results_cpd2)
            self.assertIn(value1, results_cpd2)

        except AssertionError:
            self.assertIn(valuev2, results_cpd2)
            self.assertIn(value1v2, results_cpd2)

        print ("...Testing reaction gene relationship information inputed into database")
        query = """SELECT gene_ID FROM reaction_gene WHERE
                reaction_ID='rxn2_c0' and model_ID='t1'"""
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], '(g2)')

        query = """SELECT gene_ID FROM reaction_gene WHERE
                reaction_ID='rxn5_c0' and model_ID='t2' """
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], '(g5 or g5b)')

        print ("...Testing reaction protein relationship information inputed into database")
        query = """SELECT protein_ID FROM reaction_protein WHERE
                reaction_ID='rxn5_c0' and model_ID='t2'"""
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], '(p5 or p5b)')


        query = """SELECT protein_ID FROM reaction_protein WHERE
                reaction_ID='rxn11_c0' and model_ID='t1'"""
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], '(p11)')
        os.remove(PATH+'/test.db')
if __name__ == '__main__':
    unittest.main()
