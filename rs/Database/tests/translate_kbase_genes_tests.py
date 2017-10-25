from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on gene translation'
import os
import unittest
from Database import translate_kbase_genes as tkg
from Database import query as Q
from Database import generate_database as gen_db
PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/kbasetest.db') is True:
    os.remove(PATH+'/kbasetest.db')
gen_db.Createdb(PATH+'/kbasetest.db', PATH+'/datam', False, 'bio')
DB = Q.Connector(PATH+'/kbasetest.db')

class Translate_kbase_genes_Tests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.G = tkg.TranslateKbaseGenes(DB, 'tests/metacyc_data/KbaseAllGenes_test.txt', OUTPUTPATH=PATH)
    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_translate_kbase_genes(self):
        self.assertIn('kb|g.215461.CDS.481', self.G.translation_dict.keys())
        self.assertEqual(self.G.translation_dict['kb|g.215461.CDS.481'],
                         'Ubiquinol-cytochrome C reductase iron-sulfur subunit (EC 1.10.2.2)')

        os.remove(PATH+'/kbasetest.db')
        os.remove(PATH+'/GeneTranslations.txt')

if __name__ == '__main__':
    unittest.main()
    