from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Modules for building FBA model from metabolic database'

from cobra import Reaction, Metabolite
def verbose_print(verbose, line):
    if verbose:
        print(line)


def load_reactions(model, model_id, rxns, media_constraints, compounds_dict, DB, verbose):
    '''Inserts necessary reactions in to CobraPy FBA model'''
    verbose_print(verbose, 'STATUS: Getting organism reactions ... ')
    for rxn in rxns:
        if rxn == 'biomass':
            rxn = rxn+'_'+model_id
        reaction = Reaction(rxn)
        reaction.name = DB.get_reaction_name(rxn)
        revers = DB.is_reversible_all(rxn)
        if len(media_constraints) > 0:
            if rxn in media_constraints:
                reaction.lower_bound = float(media_constraints[rxn][0])
                reaction.upper_bound = float(media_constraints[rxn][1])
            elif rxn.startswith('EX'):
                reaction.lower_bound = float(0)
                reaction.upper_bound = float(1000)
            else:
                if revers == 'true':
                    reaction.lower_bound = float(-1000)
                    reaction.upper_bound = float(1000)
                else:
                    reaction.lower_bound = float(0)
                    reaction.upper_bound = float(1000)
        else:
            if revers == 'true':
                reaction.lower_bound = float(-1000)
                reaction.upper_bound = float(1000)
            else:
                reaction.lower_bound = float(0)
                reaction.upper_bound = float(1000)
        reactants = DB.get_reactants(rxn)
        products = DB.get_products(rxn)
        met_dict = {}
        if reactants:
            for react in reactants:
                s = DB.get_stoichiometry(rxn, react, 0)
                met_dict[model.metabolites.get_by_id(compounds_dict[react])] = -float(s[0])
        if products:
            for prod in products:
                s = DB.get_stoichiometry(rxn, prod, 1)
                met_dict[model.metabolites.get_by_id(compounds_dict[prod])] = float(s[0])
        reaction.add_metabolites(met_dict)
        if DB.get_genes(rxn, model_id) is not None:
            reaction.gene_reaction_rule = DB.get_genes(rxn, model_id)
        model.add_reaction(reaction)
    return model

def load_compounds(model, metabolites, DB, verbose):
    '''Inserts necessary compounds into CobraPy FBA model'''
    verbose_print(verbose, 'STATUS: Getting organism compounds ... ')
    compounds = {}
    for c, metabolite in enumerate(metabolites):
        compounds[metabolite] = 'cpd'+str(c)
        met = Metabolite('cpd'+str(c), name='cpd'+str(c),
                         compartment=DB.get_compound_compartment(metabolite))
        model.add_metabolites(met)

    return(model, compounds)
def load_media(media_constraints, media_file):
    '''Loads specified media'''
    filename = open(media_file)
    line = filename.readline()
    if line.startswith('compounds'):
        pass
    for line in filename:
        larray = line.strip('\n').split('\t')
        media_constraints['EX_'+larray[0]+'_e0'] = [larray[2], larray[3]]
    return media_constraints
