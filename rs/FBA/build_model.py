from __future__ import print_function
__author__ = 'Leanne Whitmore and Lucy Chian'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'build flux balance analysis simulation using cobrapy'

import os
from cobra import Model
#from cobra.flux_analysis.loopless import construct_loopless_model
from FBA.generating_model import generate_model_components as gmc
PATH = os.path.dirname(os.path.abspath(__file__))

def verbose_print(verbose, line):
    if verbose:
        print(line)

class BuildModel(object):
    """
    Constructs a CobraPy FBA metabolic model for aspecific organism and runs FBA on new model
    """
    def __init__(self, target_organism_ID, inmets, inrxns, db, verbose, media=None):
        '''Initialize class'''
        self.verbose = verbose
        verbose_print(self.verbose, 'STATUS: Building FBA model ... ')
        self.target_org = target_organism_ID
        self.model = Model(target_organism_ID)
        self.rxns = inrxns
        self.metabolites = list(set(inmets))
        self.DB = db
        self.media_constraints = {}
        if media != 'Complete' and media:
            verbose_print(self.verbose, 'STATUS: loading glucose media')
            self.media_constraints = gmc.load_media(self.media_constraints, PATH+'/media/{}.tsv'.format(media))
        elif media == 'Complete' or not media:
            pass
        self.build_model()
        self.run_flux()

    def build_model(self):
        '''
        Adds compounds and rxns to cobrapy FBA model
        '''
        self.model, self.compounds_dict = gmc.load_compounds(self.model, self.metabolites, self.DB, self.verbose)
        self.model = gmc.load_reactions(self.model, self.target_org, self.rxns,
                                        self.media_constraints, self.compounds_dict,
                                        self.DB, self.verbose)

    def set_objective_function(self, biomass_rxn):
        biomass = self.model.reactions.get_by_id(biomass_rxn)
        objective_dict = {biomass:1}
        self.model.objective = objective_dict

    def run_flux(self):
        '''
        Runs FBA on model
        '''
        try:
            self.set_objective_function('biomass0_'+self.target_org)
        except KeyError:
            try:
                self.set_objective_function('bio10')
            except KeyError:
                try:
                    self.set_objective_function('bio10_'+self.target_org)
                except KeyError:   
                    try:
                        self.set_objective_function('bio1')
                    except KeyError:
                        print ('WARNING: No biomass rxn')
        #print 'STATUS: Construct loopless model'
        #self.solution1=construct_loopless_model(self.model).optimize()
        #print self.solution1
        verbose_print(self.verbose, 'STATUS: Finished optimizing first model')
        self.solution = self.model.optimize()
