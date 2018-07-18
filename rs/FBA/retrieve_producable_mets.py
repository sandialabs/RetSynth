from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'get metabolites that can be produced'

import os
from cobra import Model, Reaction
from tqdm import tqdm
from FBA.generating_model import generate_model_components as gmc
PATH = os.path.dirname(os.path.abspath(__file__))

def RetrieveActiveRxnsCompounds(target_organism_ID, inmets, inrxns, DB, output_queue, verbose, media=None):
    """Identifies active compounds in an organism"""
    media_constraints = {}
    active_metabolism = {}
    if media:
        media_constraints = gmc.load_media(media_constraints, media)
    else:
        print ('STATUS: loading glucose media')
        media_constraints = gmc.load_media(media_constraints, PATH+'/KBaseMedia_Carbon-D-Glucose_MediaCompounds.tsv')

    inmets = list(set(inmets))
    removecpds = []
    removerxns = {}
    finalremoverxns = []
    model = Model(target_organism_ID)
    model, compound_dict = gmc.load_compounds(model, inmets, DB, verbose)
    model = gmc.load_reactions(model, target_organism_ID, inrxns,
                               media_constraints, compound_dict, DB, verbose)

    print ('STATUS: Determine what compounds in network are active (can be produced)')
    for cpd in tqdm(inmets):
        model = add_sink_reaction(cpd, compound_dict, model)
        sink = model.reactions.get_by_id('Sink_'+compound_dict[cpd])
        objective_dict = {sink:1}
        model.objective = objective_dict
        solution = model.optimize()
        removerxns, removecpds = parse_solution(solution, cpd, removerxns, removecpds)

    for rxn in removerxns:
        if len(removecpds) == len(removerxns[rxn]):
            finalremoverxns.append(rxn)

    inmets = list(set(inmets)-set(removecpds))
    inrxns = list(set(inrxns)-set(finalremoverxns))
    active_metabolism[target_organism_ID] = []
    active_metabolism[target_organism_ID].append(inmets)
    active_metabolism[target_organism_ID].append(inrxns)
    if output_queue:
        output_queue.put(active_metabolism)
    else:
        return active_metabolism

def parse_solution(solution, cpd, removerxns, removecpds):
    '''
    Determines if compound can be produced
    '''
    if round(solution.f, 2) == 0:
        removecpds.append(cpd)
        for rxn in solution.x_dict.keys():
            if round(solution.x_dict[rxn], 2) == 0:
                removerxns.setdefault(rxn, []).append(cpd)
    return (removerxns, removecpds)

def add_sink_reaction(sink_cpd, compound_dict, model):
    '''
    Adds reaction for optimizing production of a specified target compound,
    used to determine if compound can be produced or not on a provided media
    '''
    reaction = Reaction('Sink_'+compound_dict[sink_cpd])
    reaction.name = 'Sink_'+compound_dict[sink_cpd]
    reaction.lower_bound = float(0)
    reaction.upper_bound = float(1000)
    met_dict_sink = {}
    met_dict_sink[model.metabolites.get_by_id(compound_dict[sink_cpd])] = float(-1)
    reaction.add_metabolites(met_dict_sink)
    model.add_reaction(reaction)
    return model
