from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Insert reactions necessary for producing target \
                    compound and run flux balance analysis'

from copy import deepcopy
from sys import exit
from cobra import Reaction, Metabolite
from tqdm import tqdm
#from cobra.flux_analysis.loopless import construct_loopless_model

def verbose_print(verbose, line):
    if verbose:
        print(line)

class OptimizeTarget(object):
    """
    Adds compounds and reactions needed to produce target compound
    into previously constructed CobraPy FBA model and runs FBA
    """
    def add_external_metabolites(self):
        '''
        Adds compounds that need to be present to produce target compound
        '''
        for count, metabolite in enumerate(self.external_metabolites):
            if metabolite in self.compounds_dict.keys():
                pass
            else:
                self.compounds_dict[metabolite] = 'E'+str(count)
                met = Metabolite('E'+str(count), name=self.DB.get_compound_name(metabolite),
                                 compartment=self.DB.get_compound_compartment(metabolite))
                self.model.add_metabolites(met)

    def add_external_reactions(self):
        '''
        Adds reactions that need to be present to produce target compound
        '''
        reactants = []
        products = []
        for rxn in self.exrxns:
            if rxn in self.inrxns:
                pass
            else:
                reaction = Reaction(rxn)
                reaction.name = self.DB.get_reaction_name(rxn)
                met_dict = {}
                if self.exrxns[rxn]['direction'] == 'reverse':
                    for react in self.exrxns[rxn]['reactants']:
                        products.append(react)
                        s = self.DB.get_stoichiometry(rxn, react, 0)
                        met_dict[self.model.metabolites.get_by_id(
                            self.compounds_dict[react])] = float(s[0])
                    for prod in self.exrxns[rxn]['products']:
                        reactants.append(prod)
                        s = self.DB.get_stoichiometry(rxn, prod, 1)
                        met_dict[self.model.metabolites.get_by_id(
                            self.compounds_dict[prod])] = -float(s[0])
                else:
                    for react in self.exrxns[rxn]['reactants']:
                        reactants.append(react)
                        s = self.DB.get_stoichiometry(rxn, react, 0)
                        met_dict[self.model.metabolites.get_by_id(
                        	   self.compounds_dict[react])] = float(-s[0])
                    for prod in self.exrxns[rxn]['products']:
                        products.append(prod)
                        s = self.DB.get_stoichiometry(rxn, prod, 1)
                        met_dict[self.model.metabolites.get_by_id(
                        	   self.compounds_dict[prod])] = float(s[0])
                reaction.add_metabolites(met_dict)
                try:
                    self.model.add_reaction(reaction)
                except ValueError:
                    pass
        sink_comp = []
        self.sink_rxns = []

        for c in products:
            if c not in reactants and c not in sink_comp and c not in self.inmets:
                sink_comp.append(c)

        for sc in sink_comp:
            self.add_sink_reaction(sc)

        if self.target_compound not in sink_comp:
            self.add_sink_reaction(self.target_compound)

    def add_sink_reaction(self, sink_cpd):
        reaction = Reaction('Sink_'+self.compounds_dict[sink_cpd])
        self.sink_rxns.append('Sink_'+self.compounds_dict[sink_cpd])
        reaction.name = 'Sink_'+self.compounds_dict[sink_cpd]
        reaction.lower_bound = float(0)
        reaction.upper_bound = float(1000)
        met_dict_sink = {}
        met_dict_sink[self.model.metabolites.get_by_id(
        	   self.compounds_dict[sink_cpd])] = float(-1)
        reaction.add_metabolites(met_dict_sink)
        self.model.add_reaction(reaction)

    def set_objective_function(self, biomass_rxn):
        biomass = self.model.reactions.get_by_id(biomass_rxn)
        sink = self.model.reactions.get_by_id('Sink_'+self.compounds_dict[self.target_compound])
        self.objective_dict = {biomass:1, sink:1}
        self.model.objective = self.objective_dict
 
    def run_fba(self):
        '''
        Runs FBA on CobraPy FBA model with added reactions and compounds
        optimizing production of target compound and (maybe) biomass
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
                        sink = self.model.reactions.get_by_id('Sink_'+self.compounds_dict[self.target_compound])
                        self.objective_dict = {sink:1}
                        self.model.objective = self.objective_dict
        # self.fbasol=construct_loopless_model(self.model).optimize()
        self.fbasol = self.model.optimize()

    def remove_modifications(self):
        '''
        Removes external reactions and compounds from an FBA model
        '''
        [self.model.metabolites.remove(self.compounds_dict[m])
         for m in self.external_metabolites]
        [self.model.remove_reactions(r) for r in self.exrxns]
        [self.model.remove_reactions(s) for s in self.sink_rxns]

    def perform_knockouts(self):
        '''
        Removes each reaction one at a time and runs FBA to see if there is
        a difference in production of target compound
        '''
        self.KOsolutions = {}
        self.objrxns_KO = {}
        self.essentialrxns = []
        self.KOfbasol = self.model.optimize()
        print ("STATUS: Performing knockouts")
        for rxn in tqdm(self.model.reactions):
            if rxn.id == 'Sink_'+self.target_compound or rxn.id == 'biomass0_'+self.target_org:
                pass
            else:
                org_lb = deepcopy(rxn.lower_bound)
                org_up = deepcopy(rxn.upper_bound)
                rxn.lower_bound = float(0)
                rxn.upper_bound = float(0)
                self.KOfbasol = self.model.optimize()
                count = 0
                for rxn_obj in self.objective_dict:
                    if (self.fbasol.fluxes[rxn_obj.id] < self.KOfbasol.fluxes[rxn_obj.id] or
                            self.fbasol.fluxes[rxn_obj.id] > self.KOfbasol.fluxes[rxn_obj.id]):
                        count += 1
                        self.objrxns_KO[rxn.id] = {}
                        self.objrxns_KO[rxn.id][rxn_obj.id] = '\t'.join([str(self.fbasol.fluxes[rxn_obj.id]),
                                                                         str(self.KOfbasol.fluxes[rxn_obj.id])])
                if count != 0:
                    self.KOsolutions[rxn.id] = self.KOfbasol
                if (round(self.fbasol.fluxes['Sink_'+self.compounds_dict[self.target_compound]], 2) >
                        round(self.KOfbasol.fluxes['Sink_'+self.compounds_dict[self.target_compound]], 2)):
                    self.essentialrxns.append(rxn.id)
                rxn.lower_bound = org_lb
                rxn.upper_bound = org_up

    def __init__(self, target_compound_ID, target_organism_ID, model, temp_rxns,
                 temp_exmets, compounds_dict, inmets, inrxns, db, verbose, KO=True, remove=True):
        '''Initialize class'''
        verbose_print(verbose, 'STATUS: Running FBA to optimize target {}...'.format(target_compound_ID))
        self.verbose = verbose
        self.inmets = inmets
        self.inrxns = inrxns
        self.target_org = target_organism_ID
        self.target_compound = target_compound_ID
        self.model = model
        self.exrxns = {}
        self.external_metabolites = []
        self.DB = db
        self.compounds_dict = compounds_dict
        for key, os_dict in temp_rxns.iteritems():
            for r in os_dict:
                if r not in self.exrxns:
                    self.exrxns[r] = os_dict[r]
        for key, exmets in temp_exmets.iteritems():
            for c in exmets:
                if c not in self.external_metabolites:
                    self.external_metabolites.append(c)
        if self.exrxns:
            self.add_external_metabolites()
            self.add_external_reactions()
            self.run_fba()
            if KO:
                self.perform_knockouts()
            if remove:
                self.remove_modifications()
        else:
            print ('ERROR: No external reactions in solution, likely issue else where in code')
            exit()
