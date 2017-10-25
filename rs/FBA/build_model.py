from __future__ import print_function
__author__ = 'Leanne Whitmore and Lucy Chian'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'build flux balance analysis simulation using cobrapy'

from cobra import Model
from cobra.flux_analysis.loopless import construct_loopless_model
from FBA.generating_model import generate_model_components as gmc

class BuildModel(object):
    """
    Constructs a CobraPy FBA metabolic model for aspecific organism and runs FBA on new model
    """
    def __init__(self, target_organism_ID, inmets, inrxns, db, media=None):
        '''Initialize class'''
        print ('STATUS: Building FBA model ... ')
        self.target_org = target_organism_ID
        self.model = Model(target_organism_ID)
        self.rxns = inrxns
        self.metabolites = inmets
        self.DB = db
        self.media_constraints = {}
        if media:
            self.media_constraints = gmc.load_media(self.media_constraints, media)
        self.build_model()
        self.run_flux()

    def build_model(self):
        '''
        Adds compounds and rxns to cobrapy FBA model
        '''
        self.model, self.compounds_dict = gmc.load_compounds(self.model, self.metabolites, self.DB)
        self.model = gmc.load_reactions(self.model, self.target_org, self.rxns,
                                        self.media_constraints, self.compounds_dict, self.DB)

    def run_flux(self):
        '''
        Runs FBA on model
        '''
        try:
            biomass = self.model.reactions.get_by_id('biomass0_'+self.target_org)
            objective_dict = {biomass:1}
            self.model.objective = objective_dict

        except KeyError:
            print ('WARNING: No biomass rxn')
        #print 'STATUS: Construct loopless model'
        #self.solution1=construct_loopless_model(self.model).optimize()
        #print self.solution1
        print ('STATUS: Finished optimizing first model')
        self.solution = self.model.optimize()
