from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests construction of an FBA model from the Database'

import os
import re
import unittest
import cobra
from FBA import build_model as bm
from Database import query as Q
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/FBA/tests', '', PATH)

if os.path.isfile(PATH+'/test.db') is True:
    os.remove(PATH+'/test.db')

'''Load database'''
init_db.Createdb(PATH+'/test.db', False)
bkdb.BuildKbase(PATH+'/data', PPATH+'/Database/KbasetoKEGGCPD.txt',
                PPATH+'/Database/KbasetoKEGGRXN.txt', False,
                PATH+'/test.db', 'bio')
DB = Q.Connector(PATH+'/test.db')
inmets = DB.get_compounds_in_model('t2')
inrxns = DB.get_reactions_in_model('t2')

'''Load model directly from sbml file'''
model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
os.remove(PATH+'/test.db')

class BuildFBATests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.fba = bm.BuildModel('t2', inmets, inrxns, DB)

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
        del self.fba

    def test_num_metabolites(self):
        """Test construction of CobraPy FBA model, specifically adding metabolites"""
        '''Test number of metabolites added to model'''
        print ("Testing construction of FBA model, addition of metabolites")
        self.assertEqual(len(self.fba.model.metabolites), len(inmets))
        '''Test reactions compound is in'''
        print ("...Test reactions compound is in")
        fbareactions = self.fba.model.metabolites.cpd0.reactions
        testreactions = model.metabolites.cpdA_c0.reactions
        self.assertEqual(len(fbareactions), len(testreactions))
        test_list = []
        for reaction in fbareactions:
            test_list.append(reaction.id)
        for reaction in testreactions:
            test_list.append(reaction.id)
        self.assertEqual(len(set(test_list)), len(fbareactions))
        '''Test compartments of compound'''
        print ("...Testing compartments of compound")
        self.assertEqual(self.fba.model.metabolites.cpd0.compartment,
                         model.metabolites.cpdA_c0.compartment)

        '''Test name of compound'''
        print ("...Testing name of compound")
        self.assertEqual(self.fba.model.metabolites.cpd0.name,
                         model.metabolites.cpdA_c0.name)

    def test_reactions(self):
        """Test construction of CobraPy FBA model, specifically adding reactions"""
        '''Test number of reactions added to model '''
        print ("Testing construction of FBA model, addition of reactions")
        self.assertEqual(len(self.fba.model.reactions), len(inrxns))

        '''Test products and reactants of reaction'''
        print ('...Testing products and reactants of a reaction')
        fbareactants = self.fba.model.reactions.rxn1_c0.reactants
        fbaproducts = self.fba.model.reactions.rxn1_c0.products
        testmodelreactants = model.reactions.rxn1_c0.reactants
        testmodelproducts = model.reactions.rxn1_c0.products
        self.assertEqual(len(fbareactants), len(testmodelreactants))
        self.assertEqual(len(fbaproducts), len(testmodelproducts))

        '''Test upper bound and lower bound of reaction'''
        print ("...Test upper bound and lower bound of reaction")
        fbalowerbound = self.fba.model.reactions.rxn5_c0.lower_bound
        testlowerbound = model.reactions.rxn5_c0.lower_bound
        self.assertEqual(fbalowerbound, testlowerbound)
        fbaupperbound = self.fba.model.reactions.rxn5_c0.upper_bound
        testupperbound = model.reactions.rxn5_c0.upper_bound
        self.assertEqual(fbaupperbound, testupperbound)

        '''Test coefficents of reaction'''
        print ("...Test coefficents of reaction")
        fbareactants = self.fba.model.reactions.rxn5_c0.reactants
        fbaproducts = self.fba.model.reactions.rxn5_c0.products
        for m in fbareactants:
            fbacoef = self.fba.model.reactions.rxn5_c0.get_coefficient(m.id)
            ind = self.fba.compounds_dict.values().index(m.id)
            testcoef = model.reactions.rxn5_c0.get_coefficient(self.fba.compounds_dict.keys()[ind])
            self.assertEqual(fbacoef, testcoef)
        for m in fbaproducts:
            fbacoef = self.fba.model.reactions.rxn5_c0.get_coefficient(m.id)
            ind = self.fba.compounds_dict.values().index(m.id)
            testcoef = model.reactions.rxn5_c0.get_coefficient(self.fba.compounds_dict.keys()[ind])
            self.assertEqual(fbacoef, testcoef)

        '''Test name of reaction'''
        print ("...Test name of reaction")
        self.assertEqual(self.fba.model.reactions.rxn5_c0.name,
                         model.reactions.rxn5_c0.name)

        '''Test gene reaction rule'''
        print ("...Test gene reaction rule")
        fbagenes = self.fba.model.reactions.rxn5_c0.gene_reaction_rule
        testgenes = model.reactions.rxn5_c0.gene_reaction_rule
        self.assertEqual(fbagenes, testgenes)

        fbagenes = self.fba.model.reactions.rxn5_c0.gene_reaction_rule
        testgenes = model.reactions.rxn5_c0.gene_reaction_rule
        self.assertEqual(fbagenes, testgenes)

        fbagenes = self.fba.model.reactions.rxn7_c0.gene_reaction_rule
        testgenes = model.reactions.rxn7_c0.gene_reaction_rule
        self.assertEqual(fbagenes, '(g7a and g7b)')
        self.assertEqual(testgenes, '(g7a and g7b)')
if __name__ == '__main__':
    unittest.main()
