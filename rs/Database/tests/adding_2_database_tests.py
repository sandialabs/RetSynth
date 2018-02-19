from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes add reactions and compounds to existing database'

import os
import unittest
import sqlite3
from Database import build_kbase_db as bkdb
PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')
sqlite_database = sqlite3.connect(PATH+'/test.db')
sqlite_database.execute('''CREATE table model (ID text,file_name text)''')
sqlite_database.execute('''CREATE table compound (ID text, name text,compartment text, kegg_id)''')
sqlite_database.execute('''CREATE table compartments (ID text, name text)''')
sqlite_database.execute('''CREATE table model_compound (cpd_ID text,model_ID)''')
sqlite_database.execute('''CREATE table reaction (ID text, name text, kegg_id, type text)''')
sqlite_database.execute('''CREATE table model_reaction (reaction_ID text,
                        model_ID text, is_rev bit(1))''')
sqlite_database.execute('''CREATE table reaction_compound (reaction_ID text,
                        cpd_ID text, is_prod bit(1), stoichiometry int, filenum)''')
sqlite_database.execute('''CREATE table reaction_reversibility (reaction_ID text,
	                    is_reversible bit(1))''')
sqlite_database.execute('''CREATE table reaction_gene (reaction_ID text,
	                    model_ID text, gene_ID text)''')
sqlite_database.execute('''CREATE table reaction_protein (reaction_ID text,
	                    model_ID text, protein_ID text)''')
sqlite_database.execute('''CREATE table cluster (cluster_num text, ID text)''')

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
sqlite_database.execute('''CREATE INDEX compound_ind ON compound(ID)''')
sqlite_database.execute('''CREATE INDEX reactiongene_ind ON
                        reaction_gene(reaction_ID, model_ID)''')
sqlite_database.execute('''CREATE INDEX reactionprotein_ind ON
                         reaction_protein(reaction_ID, model_ID)''')
sqlite_database.execute('''CREATE INDEX cluster_ind ON cluster(cluster_num)''')
sqlite_database.execute('''CREATE INDEX reaction_reversibility_ind ON
                        reaction_reversibility(reaction_ID)''')
class Adding_to_current_DBTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_adding_to_db(self):
        """Tests if compounds and reactions are being added to a preexisting database correctly"""
        bkdb.BuildKbase(PATH+'/data', '../KbasetoKEGGCPD.txt', '../KbasetoKEGGRXN.txt', False, PATH+'/test.db', 'bio')

        bkdb.BuildKbase(PATH+'/data2', '../KbasetoKEGGCPD.txt', '../KbasetoKEGGRXN.txt', False, PATH+'/test.db')
        print ("...Testing reaction reversibility information getting added to a pre-existing database")
        Q = sqlite_database.execute("""SELECT is_reversible FROM reaction_reversibility
                                    WHERE reaction_ID='rxn21_c0'""")
        result = Q.fetchone()
        self.assertEqual(result[0], 'false')

        Q = sqlite_database.execute("""SELECT is_reversible FROM reaction_reversibility
                                    WHERE reaction_ID='rxn20_c0'""")
        result = Q.fetchone()
        self.assertEqual(result[0], 'true')

        print ("...Testing reaction information getting added to a pre-existing database")
        Q = sqlite_database.execute("SELECT ID FROM reaction WHERE ID='rxn21_c0'")
        result = Q.fetchone()
        self.assertEqual(result[0], 'rxn21_c0')

        Q = sqlite_database.execute("SELECT ID FROM reaction WHERE ID='rxn20_c0'")
        result = Q.fetchone()
        self.assertEqual(result[0], 'rxn20_c0')

        print ("...Testing reaction compound information getting added to a pre-existing database")
        Q = sqlite_database.execute("""SELECT reaction_ID FROM reaction_compound
                                    WHERE reaction_ID='rxn21_c0'""")
        result = Q.fetchone()
        self.assertEqual(result[0], 'rxn21_c0')

        Q = sqlite_database.execute("""SELECT count(*) FROM reaction_compound
                                     WHERE reaction_ID='rxn20_c0'""")
        result = Q.fetchone()
        self.assertEqual(result[0], 3)

        Q = sqlite_database.execute("""SELECT count(*) FROM reaction_compound
                                     WHERE reaction_ID='rxn21_c0'""")
        result = Q.fetchone()
        self.assertEqual(result[0], 2)

        print ("...Testing compound information getting added to a pre-existing database")
        Q = sqlite_database.execute("SELECT ID FROM compound WHERE ID='cpdX_c0'")
        result = Q.fetchone()
        self.assertEqual(result[0], 'cpdX_c0')

        Q = sqlite_database.execute("SELECT ID FROM compound WHERE ID='cpdW_c0'")
        result = Q.fetchone()
        self.assertEqual(result[0], 'cpdW_c0')

        Q = sqlite_database.execute("SELECT ID FROM compound WHERE ID='cpdZ_c0'")
        result = Q.fetchone()
        self.assertEqual(result[0], 'cpdZ_c0')

        Q = sqlite_database.execute("SELECT ID FROM compound WHERE ID='cpdY_c0'")
        result = Q.fetchone()
        self.assertEqual(result[0], 'cpdY_c0')

        print ("...Testing cluster_num information")
        Q = sqlite_database.execute("SELECT DISTINCT cluster_num FROM cluster")
        hits = Q.fetchall()
        uniq_clusters = [i[0] for i in hits]
        self.assertEqual(len(uniq_clusters), 3)
        os.remove(PATH+'/test.db')

if __name__ == '__main__':
    unittest.main()
