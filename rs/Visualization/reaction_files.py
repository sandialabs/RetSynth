from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Generate Reaction smiles file (can be used in ChemDraw)'
import re
import os
import httplib
import urllib2
from sys import platform
import pubchempy as pcp
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_renderer
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_renderer
elif platform == "win32" or platform == 'win64':
    raise ImportError('Cannot translate RDF file on windows machine')


class ReactionFiles(object):
    def __init__(self, output_path, DB, reactions, target,
                 target_organism, incpds, figures, cdxmlfiles=True):
        self.output_path = output_path
        self.DB = DB
        self.incpds = incpds
        self.incpds_set = set(incpds)
        self.reactions = reactions
        self.target = target
        self.target_organism = target_organism
        self.figures = figures
        self.target_organism_name = self.DB.get_organism_name(target_organism)
        self.IN = indigo.Indigo()
        self.IR = indigo_renderer.IndigoRenderer(self.IN)
        self.cdxmlfiles = cdxmlfiles
        self.generate_reaction_SMILES()

    def get_SMILES(self, compounds, smiles_cpds):
        '''Get smiles for a compounds'''
        for cpd in compounds:
            if cpd.startswith('InChI'):
                cpd = re.sub('_\w+\d+$', '', cpd)
                mol = self.IN.loadMolecule(cpd)
                SMILE = mol.canonicalSmiles()
                smiles_cpds.append(SMILE)
            else:
                cpd_name = self.DB.get_compound_name(cpd)
                cpd_name = re.sub('_\w+\d+$', '', cpd_name)
                cpd_name = re.sub('_', ' ', cpd_name)
                name_cpds = self.get_smiles_from_name(cpd_name, smiles_cpds)
                if not name_cpds:
                    cpd_name = re.sub(' ', '-', cpd_name)
                    name_cpds = self.get_smiles_from_name(cpd_name, smiles_cpds)
                    if not name_cpds:
                        cpd_name_p = cpd_name+'+'
                        name_cpds = self.get_smiles_from_name(cpd_name_p, smiles_cpds)
                        if not name_cpds:
                            cpd_name_m = cpd_name+'-'
                            name_cpds = self.get_smiles_from_name(cpd_name_m, smiles_cpds)
                if name_cpds:
                    name_cpds = smiles_cpds
                else:
                    print ('WARNING: No SMILE was identifed for {} - {}'.format(cpd, cpd_name))
        return smiles_cpds

    def get_smiles_from_name(self, cpd_name, smiles_cpds):
        '''Get smiles for a compound from its' name'''
        try:
            pcp_cpds = pcp.get_compounds(cpd_name, 'name')
            if pcp_cpds:
                pcp_cpd = pcp_cpds[0]
                SMILE = pcp_cpd.canonical_smiles
                smiles_cpds.append(SMILE)
            return smiles_cpds
        except (pcp.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError, ValueError):
            return None

    def generate_output_folders(self, folder):
        '''Check to see if output folders are present
            if not set up new output folder'''
        try:
            os.mkdir(folder)
        except OSError:
            pass

    def order_of_paths(self):
        ordered_paths = {}
        for count_pathway, os_dict in self.reactions.iteritems():
            ordered_paths[count_pathway] = []
            counter = 0
            store_reactants = []
            rxn_set = set()
            for rxn in os_dict:
                if os_dict[rxn]['direction'] == 'forward':
                    if self.target in os_dict[rxn]['products']:
                        counter = 1
                        ordered_paths[count_pathway] = {}
                        ordered_paths[count_pathway][counter] = rxn
                        rxn_set.add(rxn)
                        store_reactants = os_dict[rxn]['reactants'].keys()
                else:
                    if self.target in os_dict[rxn]['reactants']:
                        counter = 1
                        ordered_paths[count_pathway] = {}
                        ordered_paths[count_pathway][counter] = rxn
                        rxn_set.add(rxn)
                        store_reactants = os_dict[rxn]['products'].keys()

            while len(ordered_paths[count_pathway]) != len(os_dict):
                for rxn in os_dict:
                    if rxn not in rxn_set:
                        if os_dict[rxn]['direction'] == 'forward':
                            products_set = set(os_dict[rxn]['products'])
                            intersect = products_set.intersection(self.incpds_set)
                            for i in intersect:
                                products_set.remove(i)
                            if products_set:
                                overlap = products_set.intersection(store_reactants)
                                if overlap:
                                    counter += 1
                                    ordered_paths[count_pathway][counter] = rxn
                                    rxn_set.add(rxn)
                                    store_reactants.extend(os_dict[rxn]['reactants'].keys())

                        else:
                            reactants_set = set(os_dict[rxn]['reactants'])
                            intersect = reactants_set.intersection(self.incpds_set)
                            for i in intersect:
                                reactants_set.remove(i)
                            if reactants_set:
                                overlap = reactants_set.intersection(store_reactants)
                                if overlap:
                                    counter += 1
                                    ordered_paths[count_pathway][counter] = rxn
                                    rxn_set.add(rxn)
                                    store_reactants.extend(os_dict[rxn]['products'].keys())
        return ordered_paths
    def generate_reaction_SMILES(self):
        '''Produce files that have reaction smiles for each reaction in a pathway'''
        ordered_paths = self.order_of_paths()
        self.ordered_paths = ordered_paths
        if self.figures:
            target_reformat = re.sub('/', '_', self.target)
            self.generate_output_folders(self.output_path+'/solution_smiles')
            self.generate_output_folders(self.output_path+'/solution_cdxml')
            self.generate_output_folders(self.output_path+'/solution_cdxml/'+
                                         target_reformat+'_solutions')
            self.generate_output_folders(self.output_path+'/solution_figures')
            self.generate_output_folders(self.output_path+'/solution_figures/'+
                                         target_reformat+'_solutions')
            for count_pathway, os_dict in self.reactions.iteritems():
                smiles_reactants = {}
                smiles_products = {}

                with open(self.output_path+'/solution_smiles/reaction_smile_'+
                          target_reformat+'_solution_'+str(count_pathway)+
                          '.smi', 'w') as fout:
                    array = self.IN.createArray()
                    for counter in reversed(ordered_paths[count_pathway].keys()):
                        rxn = ordered_paths[count_pathway][counter]
                        smiles_reactants[rxn] = []
                        smiles_products[rxn] = []
                        if os_dict[rxn]['direction'] == 'forward':
                            smiles_reactants[rxn] = self.get_SMILES(os_dict[rxn]['reactants'],
                                                                    smiles_reactants[rxn])
                            smiles_products[rxn] = self.get_SMILES(os_dict[rxn]['products'],
                                                                   smiles_products[rxn])
                        else:
                            smiles_reactants[rxn] = self.get_SMILES(os_dict[rxn]['products'],
                                                                    smiles_reactants[rxn])
                            smiles_products[rxn] = self.get_SMILES(os_dict[rxn]['reactants'],
                                                                   smiles_products[rxn])
                    for count_rxn in reversed(ordered_paths[count_pathway].keys()):
                        key = ordered_paths[count_pathway][count_rxn]
                        fout.write('.'.join(smiles_reactants[key])+'>>'+
                                   '.'.join(smiles_products[key])+'\n')
                        rxn = self.IN.loadReaction('.'.join(smiles_reactants[key])+'>>'+
                                                   '.'.join(smiles_products[key]))
                        self.IN.setOption('render-grid-title-property', str(key)+
                                          ' '+str(self.DB.get_reaction_name(key)))
                        array.arrayAdd(rxn)
                        self.IR.renderToFile(rxn, self.output_path+'/solution_cdxml/'+target_reformat+
                                             '_solutions/'+'Solution_'+str(count_pathway)+'_rxn_'+
                                             str(count_rxn)+'_'+self.target_organism_name+'.cdxml')
                    self.IN.setOption('render-output-format', 'png')
                    self.IN.setOption("render-bond-length", "50")
                    self.IN.setOption("render-margins", "80, 30, 10, 30")
                    self.IN.setOption("render-grid-margins", "5, 5, 5, 5")
                    self.IR.renderGridToFile(array, None, len(smiles_reactants.keys()),
                                             self.output_path+'/solution_figures/'+target_reformat+
                                             '_solutions/'+'Solution_'+str(count_pathway)+'_'+
                                             self.target_organism_name+'.png')
