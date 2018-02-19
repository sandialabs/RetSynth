from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests RDF reader'
import sqlite3
import os
import unittest
from shutil import copyfile
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import build_SPRESI_db as bspresidb

PATH = os.path.dirname(os.path.abspath(__file__))

if os.path.isfile(PATH+'/testRDF.db') is True:
    os.remove(PATH+'/testRDF.db')

if os.path.isfile(PATH+'/testinchiRDF.db') is True:
    os.remove(PATH+'/testinchiRDF.db')

init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data5', '../KbasetoKEGGCPD.txt',
                '../KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
copyfile(PATH+'/test.db', PATH+'/testRDF.db')
os.remove(PATH+'/test.db')
DB = Q.Connector(PATH+'/testRDF.db')
compartmentID = DB.get_compartment('cytosol')
compartmentID = compartmentID[0]
bspresidb.RDF_Reader(PATH+'/data6/', PATH+'/testRDF.db', 'chem',compartmentID, 1, temp_option=False, pressure_option=False,
                          yield_option=False, time_option=False, catalyst_option=False, solvent_option=False)


init_db.Createdb(PATH+'/testinchi.db', False)
bkdb.BuildKbase(PATH+'/data5', '../KbasetoKEGGCPD.txt',
                '../KbasetoKEGGRXN.txt', False,
                PATH+'/testinchi.db', 'bio')
copyfile(PATH+'/testinchi.db', PATH+'/testinchiRDF.db')
os.remove(PATH+'/testinchi.db')
DBinchi = Q.Connector(PATH+'/testinchiRDF.db')
compartmentID = DBinchi.get_compartment('cytosol')
compartmentID = compartmentID[0]
bspresidb.RDF_Reader(PATH+'/data6/', PATH+'/testinchiRDF.db', 'chem', compartmentID, 1, temp_option=False, pressure_option=False,
                          yield_option=False, time_option=False, catalyst_option=False, solvent_option=False)



class RDFileReaderTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def tests_RDFileReader(self):
        '''Tests RDFReader'''
        print ("Tests RDFReader")
        cnx = sqlite3.connect(PATH+'/testRDF.db')

        '''Reaction check'''
        Q = cnx.execute("SELECT reaction_ID FROM model_reaction WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 9)
        self.assertIn(('rxn4073836_s',), result)
        self.assertIn(('rxn4073837_s',), result)
        self.assertIn(('rxn4073869_s',), result)
        self.assertIn(('rxn40738361_s',), result)
        self.assertIn(('rxn40738372_s',), result)
        self.assertIn(('rxn40738693_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4073836_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073836_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4073837_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073837_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4073869_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073869_s',), result)

        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn40738361_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738361_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn40738372_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738372_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn40738693_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738693_s',), result)

        '''Reversibility check'''
        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4073836_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073836_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4073837_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073837_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4073869_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073869_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn40738361_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738361_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn40738372_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738372_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn40738693_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738693_s', 'false'), result)

        '''Compound check'''
        Q = cnx.execute("SELECT cpd_ID FROM model_compound WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 14)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',),
                         result)
        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C11H10O/c1-12-11-7-6-9-4-2-3-5-10(9)8-11/h2-8H,1H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C11H10O/c1-12-11-7-6-9-4-2-3-5-10(9)8-11/h2-8H,1H3_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C6H7N/c7-6-4-2-1-3-5-6/h1-5H,7H2_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C6H7N/c7-6-4-2-1-3-5-6/h1-5H,7H2_c0',),
                         result)

        '''reaction_compound check'''
        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn4073836_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073836_s', 'InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0', 0, 1,
                       0), result)
        self.assertIn(('rxn4073836_s', 'InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0', 1, 1, 0),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn4073869_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073869_s', 'InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0', 0, 1, 0),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn40738361_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738361_s', 'InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0', 0, 1,
                       0), result)
        self.assertIn(('rxn40738361_s', 'InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0', 1, 1, 0),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', 'InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0', 0, 1, 0),
                      result)

        '''Catalysts check'''
        Q = cnx.execute("SELECT * FROM reaction_catalysts WHERE reaction_ID = ?",
                        ('rxn4073836_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073836_s', 'InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)_c0', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_catalysts WHERE reaction_ID = ?",
                        ('rxn40738361_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738361_s', 'InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)_c0', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_catalysts")
        result = Q.fetchall()
        self.assertEqual(len(result), 12)

        Q = cnx.execute("SELECT * FROM reaction_solvents")
        result = Q.fetchall()
        self.assertEqual(len(result), 4)

        '''Solvents check'''
        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?", ('rxn4073869_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073869_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3_c0', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?",
                        ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3_c0', 'None'), result)
        
        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?",
                        ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3_c0', 'None'), result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', '120 degree', 'None', '12 h', '93.0-93.0',
                       'JOURNAL ARTICLE Tobisu Mamoru; Shimasaki Toshiaki; Chatani Naoto; Angew. Chem., Int. Ed. +Engl., 2008, Vol. 47, P. 4866-4869'), result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn46000844_s',))
        result = Q.fetchall()
        self.assertIn(('rxn46000844_s', '25 degree,80 degree,None, 600 degree', 'None,None,None,None',
                       '1 h, 30 min,None,None', '44.0-44.0',
                       'JOURNAL ARTICLE Begue Didier; Dargelos Alain; Berstermann Hans M.; Netsch Klaus P.; Bedna+rek Pawel; Wentrup Curt; J. Org. Chem., 2014, Vol. 79, P. 1247-1253'), result)
    
        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn2246995_s',))
        result = Q.fetchall()
        self.assertIn(('rxn2246995_s', 'None', 'None', 'None', 'None',
                       'PATENT CATALYST AND PROCESS FOR THE SELECTIVE DIMERIZATION OF PROPYLENE TO METHY+L-1-PENTENE, STEVENS JAMES C.; FORDYCE WILLIAM A., 1991, Patent number-4855523, US, Patent Class-4 C 07 C 2/10, Patent Owner-THE DOW CHEMICAL CO.'), result)

        os.remove(PATH+'/testRDF.db')


    def tests_RDFileReader_inchi(self):
        '''Tests RDFReader with inchi database'''
        print ("Tests RDFReader with inchi database")
        cnx = sqlite3.connect(PATH+'/testinchiRDF.db')

        '''Reaction check'''
        Q = cnx.execute("SELECT reaction_ID FROM model_reaction WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 9)
        self.assertIn(('rxn4073836_s',), result)
        self.assertIn(('rxn4073837_s',), result)
        self.assertIn(('rxn4073869_s',), result)
        self.assertIn(('rxn40738361_s',), result)
        self.assertIn(('rxn40738372_s',), result)
        self.assertIn(('rxn40738693_s',), result)

        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4073836_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073836_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4073837_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073837_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4073869_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073869_s',), result)

        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn40738361_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738361_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn40738372_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738372_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn40738693_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738693_s',), result)

        '''Reversibility check'''
        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4073836_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073836_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4073837_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073837_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4073869_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4073869_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn40738361_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738361_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn40738372_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738372_s', 'false'), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn40738693_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn40738693_s', 'false'), result)

        '''Compound check'''
        Q = cnx.execute("SELECT cpd_ID FROM model_compound WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 14)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C11H10O/c1-12-11-7-6-9-4-2-3-5-10(9)8-11/h2-8H,1H3_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C11H10O/c1-12-11-7-6-9-4-2-3-5-10(9)8-11/h2-8H,1H3_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('InChI=1S/C6H7N/c7-6-4-2-1-3-5-6/h1-5H,7H2_c0',))
        result = Q.fetchone()
        self.assertEqual(('InChI=1S/C6H7N/c7-6-4-2-1-3-5-6/h1-5H,7H2_c0',), result)

        '''reaction_compound check'''
        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn4073836_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073836_s', 'InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0', 0, 1,
                       0), result)
        self.assertIn(('rxn4073836_s', 'InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0', 1, 1, 0),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn4073869_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073869_s', 'InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0', 0, 1, 0),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn40738361_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738361_s', 'InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3_c0', 0, 1,
                       0), result)
        self.assertIn(('rxn40738361_s', 'InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3_c0', 1, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', 'InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3_c0', 0, 1, 0), result)

        '''Catalysts check'''
        Q = cnx.execute("SELECT * FROM reaction_catalysts WHERE reaction_ID = ?", ('rxn4073836_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073836_s', 'InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)_c0', 'None'),
                      result)
    
        Q = cnx.execute("SELECT * FROM reaction_catalysts WHERE reaction_ID = ?",
                        ('rxn40738361_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738361_s', 'InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)_c0', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_catalysts")
        result = Q.fetchall()
        self.assertEqual(len(result), 12)

        Q = cnx.execute("SELECT * FROM reaction_solvents")
        result = Q.fetchall()
        self.assertEqual(len(result), 4)

        '''Solvents check'''
        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?", ('rxn4073869_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4073869_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3_c0', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?",
                        ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3_c0', 'None'),
                      result)
        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn40738693_s',))
        result = Q.fetchall()
        self.assertIn(('rxn40738693_s', '120 degree', 'None', '12 h', '93.0-93.0',
                       'JOURNAL ARTICLE Tobisu Mamoru; Shimasaki Toshiaki; Chatani Naoto; Angew. Chem., Int. Ed. +Engl., 2008, Vol. 47, P. 4866-4869'), result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn46000844_s',))
        result = Q.fetchall()
        self.assertIn(('rxn46000844_s', '25 degree,80 degree,None, 600 degree', 'None,None,None,None',
                       '1 h, 30 min,None,None', '44.0-44.0',
                       'JOURNAL ARTICLE Begue Didier; Dargelos Alain; Berstermann Hans M.; Netsch Klaus P.; Bedna+rek Pawel; Wentrup Curt; J. Org. Chem., 2014, Vol. 79, P. 1247-1253'), result)
        
        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn2246995_s',))
        result = Q.fetchall()
        self.assertIn(('rxn2246995_s', 'None', 'None', 'None', 'None',
                       'PATENT CATALYST AND PROCESS FOR THE SELECTIVE DIMERIZATION OF PROPYLENE TO METHY+L-1-PENTENE, STEVENS JAMES C.; FORDYCE WILLIAM A., 1991, Patent number-4855523, US, Patent Class-4 C 07 C 2/10, Patent Owner-THE DOW CHEMICAL CO.'), result)
        os.remove(PATH+'/testinchiRDF.db')

if __name__ == '__main__':
    unittest.main()
