from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that generate database on test set of metabolic networks'

import os
import sqlite3
import unittest
from Database import initialize_database as init_db
'''Build test database'''
PATH = os.path.dirname(os.path.abspath(__file__))

class Generate_database(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_generate_db(self):
        """Tests to see if database is filling model table with correct number of models"""
        print ("Tests to see if database is generating tables and indexes")
        if os.path.isfile(PATH+'/test.db') is True:
            os.remove(PATH+'/test.db')
        init_db.Createdb(PATH+'/test.db', False)

        query = """select count(*) from sqlite_master where type='table' and name='reaction'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='reaction_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reaction_reversibility_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table' and name='model'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='model_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table' and name='compound'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='compound_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='model_compound'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelcompound_ind1'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelcompound_ind2'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='model_reaction'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelreaction_ind1'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelreaction_ind2'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_compound'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactioncompound_ind1'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactioncompound_ind2'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_reversibility'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_gene'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactiongene_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_protein'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactionprotein_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table' and name='cluster'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='cluster_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)
        #os.remove(PATH+'/test.db')

    def test_generate_db_index(self):
        print ("""Tests to see if database is generating tables and indexes if
                      inchi option is specified""")
        if os.path.isfile(PATH+'/test.db') is True:
            os.remove(PATH+'/test.db')
        init_db.Createdb(PATH+'/test.db', True)

        query = """select count(*) from sqlite_master where type='table' and name='reaction'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='reaction_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reaction_reversibility_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table' and name='model'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='model_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table' and name='compound'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='compound_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='model_compound'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelcompound_ind1'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelcompound_ind2'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='model_reaction'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelreaction_ind1'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='modelreaction_ind2'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_compound'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactioncompound_ind1'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactioncompound_ind2'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_reversibility'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_gene'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactiongene_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='reaction_protein'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='reactionprotein_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table' and name='cluster'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index' and name='cluster_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='table'
                and name='original_db_cpdIDs'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)

        query = """select count(*) from sqlite_master where type='index'
                and name='original_db_cpdIDs_ind'"""
        cnx = sqlite3.connect(PATH+'/test.db')
        Q = cnx.execute(query,)
        hit = Q.fetchone()
        self.assertEqual(hit[0], 1)
        os.remove(PATH+'/test.db')

if __name__ == '__main__':
    unittest.main()
