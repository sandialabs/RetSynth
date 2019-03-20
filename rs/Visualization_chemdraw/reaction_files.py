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
    from indigopython130_win import indigo
    from indigopython130_win import indigo_renderer
from exclude_cpds import retrieve_promiscuous_mets
from cdxml_editor import CDXML_Editor


''' Tree Structure for Ordered Path '''
class Node(object):
    def __init__(self, root=None):
        self.root = root
        self.children = []

    def addChild(self, child):
        self.children += [child]

class Tree(object):
    def __init__(self, root=None):
        self.root = root

    def set_root(self, new_root):
        self.root = new_root

    def list_nodes(self):
        ''' Depth first search, outputs list of all nodes'''
        def traverse(root):
            output = []
            for child in root.children:
                output += traverse(child)
            output += [root.root]
            return output

        return traverse(self.root)


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
        self.ordered_paths = self.order_of_paths()
        self.promiscuous = retrieve_promiscuous_mets(DB)
    # def get_SMILES(self, compounds, smiles_cpds):
    #     '''Get smiles for a compounds'''
    #     for cpd in compounds:
    #         cpd_name = self.DB.get_compound_name(cpd)
    #         if cpd.startswith('InChI'):                
    #             cpd = re.sub('_\w+\d+$', '', cpd)
    #             mol = self.IN.loadMolecule(cpd)               
    #             SMILE = mol.canonicalSmiles()
    #             smiles_cpds.append(SMILE)
    #         else:
    #             cpd_name = re.sub('_\w+\d+$', '', cpd_name)
    #             cpd_name = re.sub('_', ' ', cpd_name)
    #             name_cpds = self.get_smiles_from_name(cpd_name, smiles_cpds)
    #             if not name_cpds:
    #                 cpd_name = re.sub(' ', '-', cpd_name)
    #                 name_cpds = self.get_smiles_from_name(cpd_name, smiles_cpds)
    #                 if not name_cpds:
    #                     cpd_name_p = cpd_name+'+'
    #                     name_cpds = self.get_smiles_from_name(cpd_name_p, smiles_cpds)
    #                     if not name_cpds:
    #                         cpd_name_m = cpd_name+'-'
    #                         name_cpds = self.get_smiles_from_name(cpd_name_m, smiles_cpds)
    #             if name_cpds:
    #                 name_cpds = smiles_cpds
    #             else:
    #                 print ('WARNING: No SMILE was identifed for {} - {}'.format(cpd, cpd_name))
    #     return smiles_cpds


    def get_cdxml(self, compounds, cdxml_cpds, promiscuous_cpds):
        '''Get cdxml files for relevant compounds'''
        for cpd in compounds:
            cpd_name = compounds[cpd]           
            if cpd.startswith('InChI'):                
                cpd = re.sub('_\w+\d+$', '', cpd)
                if cpd_name == 'None':
                    cpd_name = '_'.join(cpd.split('/'))

                mol = self.IN.loadMolecule(cpd)                              
                if self.promiscuous.get(cpd_name) is None:
                    cdxml_cpds.append(cpd_name)
                    self.IN.setOption("render-comment", cpd_name)
                    self.IN.setOption("render-output-format", "cdxml")
                    self.IR.renderToFile(mol, self.output_path+'/compounds/'+cpd_name+'.cdxml')
                else:
                    promiscuous_cpds.append(cpd_name)
            else:
                if cpd_name == 'None':
                    cpd_name = cpd
                cdxml_cpds.append(cpd_name)
        return cdxml_cpds, promiscuous_cpds


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
        ''' Ordered pathway with target as tree root '''
        ordered_paths = {}

        for count_pathway, os_dict in self.reactions.iteritems():
            ordered_paths[count_pathway] = Tree()
            rxn_set = set(os_dict)
            node_dict = {}

            for rxn in rxn_set:
                if os_dict[rxn]['direction'] == 'forward' and self.target in os_dict[rxn]['products']:
                    node_dict[rxn] = Node(rxn)
                    ordered_paths[count_pathway].set_root(node_dict[rxn])
                    rxn_set.remove(rxn)
                    break
                elif os_dict[rxn]['direction'] == 'reverse' and self.target in os_dict[rxn]['reactants']:
                    node_dict[rxn] = Node(rxn)
                    ordered_paths[count_pathway].set_root(node_dict[rxn])
                    rxn_set.remove(rxn)
                    break
            

            def add_all_children(node, rxn_set):
                ''' Adds all children of root and removes rxn from rxn_set '''
                if len(rxn_set) == 0:
                    return

                store_reactants = []
                if os_dict[node.root]['direction'] == 'forward':
                    store_reactants = os_dict[node.root]['reactants'].keys()
                else:
                    store_reactants = os_dict[node.root]['products'].keys()

                for sr in store_reactants:
                    for rxn in rxn_set:
                        if os_dict[rxn]['direction'] == 'forward' and sr in os_dict[rxn]['products']:
                            node_dict[rxn] = Node(rxn)
                            node.addChild(node_dict[rxn])
                            rxn_set.remove(rxn)
                            add_all_children(node_dict[rxn], rxn_set)
                            break
                        elif os_dict[rxn]['direction'] == 'reverse' and sr in os_dict[rxn]['reactants']:
                            node_dict[rxn] = Node(rxn)
                            node.addChild(node_dict[rxn])
                            rxn_set.remove(rxn)
                            add_all_children(node_dict[rxn], rxn_set)
                            break
                    if len(rxn_set) == 0:
                        break

            add_all_children(ordered_paths[count_pathway].root, rxn_set)

        return ordered_paths


    def alter_name_length(self, path_to_figure, cpdname):
        '''Shorten compound name if it is too long'''
        if len(path_to_figure) > 250:
            remove_variable = len(path_to_figure) - 250
            cpdname = cpdname[:-remove_variable]
        return cpdname


    def generate_cdxml_files(self, fba_fluxes=None, show_rxn_info=True):
        '''Produce cdxml files for each pathway'''   
        _FBA = fba_fluxes is not None

        target_reformat_orig = re.sub('/', '_', self.target)
        self.generate_output_folders(self.output_path+'/solution_figures')                        
        self.generate_output_folders(self.output_path+'/solution_figures/'+ 
            target_reformat_orig+'_solutions')
        self.generate_output_folders(self.output_path+'/compounds')

        # For each pathway
        for count_pathway, os_dict in self.reactions.iteritems():
            MAIN_CE = CDXML_Editor(output_path=self.output_path + '/solution_figures/' + target_reformat_orig +
                '_solutions/solution_' + str(count_pathway) + '.cdxml')

            path = self.ordered_paths[count_pathway]

            cdxml_reactants = {}
            cdxml_products = {}
            promiscuous_reactants = {}
            promiscuous_products = {}
            misc_products = {}

            reaction_proteins = {}
            reaction_solvents = {}
            reaction_catalysts = {}
            reaction_SPRESI_info = {}

            fba_values = {}

            def get_info(root, parent=None):
                for child in root.children:
                    get_info(child, parent=root)

                rxn = root.root
                cdxml_reactants[rxn] = []
                cdxml_products[rxn] = []
                promiscuous_reactants[rxn] = []
                promiscuous_products[rxn] = []
                misc_products[rxn] = []

                reaction_proteins[rxn] = None
                reaction_solvents[rxn] = None
                reaction_catalysts[rxn] = None
                reaction_SPRESI_info[rxn] = None


                if os_dict[rxn]['direction'] == 'forward':
                    cdxml_reactants[rxn], promiscuous_reactants[rxn] = self.get_cdxml(os_dict[rxn]['reactants'],
                            cdxml_reactants[rxn], promiscuous_reactants[rxn])
                    cdxml_products[rxn], promiscuous_products[rxn] = self.get_cdxml(os_dict[rxn]['products'],
                            cdxml_products[rxn], promiscuous_products[rxn])
                else:
                    cdxml_reactants[rxn], promiscuous_reactants[rxn] = self.get_cdxml(os_dict[rxn]['products'],
                            cdxml_reactants[rxn], promiscuous_reactants[rxn])
                    cdxml_products[rxn], promiscuous_products[rxn] = self.get_cdxml(os_dict[rxn]['reactants'],
                            cdxml_products[rxn], promiscuous_products[rxn])

                if _FBA:
                    fba_values[rxn] = float(max(0,fba_fluxes[rxn]))              

                protein = self.DB.get_proteins(rxn, self.target_organism)
                if protein != "None":
                    reaction_proteins[rxn] = protein                    
                elif show_rxn_info:
                    # RXN SOLVENT
                    solvents = self.DB.get_solvents(rxn)
                    if solvents != "None" and len(solvents) > 0:
                        reaction_solvents[rxn] = "Solvent: %s" % ','.join(solvents)
                    # RXN CATALYST
                    catalysts = self.DB.get_catalysts(rxn)
                    if catalysts != "None" and len(catalysts) > 0:
                        reaction_catalysts[rxn] = "Catalyst: %s" % ','.join(catalysts)


                    # RXN SPRESI
                    spresi_info = []
                    rxn_temperatures = self.DB.get_temperature(rxn)
                    rxn_pressures = self.DB.get_pressure(rxn)
                    rxn_times = self.DB.get_time(rxn)
                    rxn_yields = self.DB.get_yield(rxn)
                    if rxn_temperatures != "None" and len(rxn_temperatures) > 0 and not all(t=='None' for t in rxn_temperatures):
                        spresi_info += ["Temperature: %s" % ', '.join(rxn_temperatures)]
                    if rxn_pressures != "None" and len(rxn_pressures) > 0 and not all(p=='None' for p in rxn_pressures):
                        spresi_info += ["Pressure: %s" % ', '.join(rxn_pressures)]
                    if rxn_times != "None" and len(rxn_times) > 0 and not all(t=='None' for t in rxn_times):
                        spresi_info += ["Time: %s" % ', '.join(rxn_times)]
                    if rxn_yields != "None" and len(rxn_yields) > 0 and not all(y=='None' for y in rxn_yields):
                        spresi_info += ["Yield: %s" % ', '.join(rxn_yields)]
                    if len(spresi_info) > 0:
                        reaction_SPRESI_info[rxn] = '\n'.join(spresi_info)
                

                if os_dict[rxn]['direction'] == 'forward':
                    cur_prods = os_dict[rxn]['products']
                else:
                    cur_prods = os_dict[rxn]['reactants']
                
                difference = []
                if parent:
                    next_rxn = parent.root 
                    if os_dict[next_rxn]['direction'] == 'forward':
                        next_reacts = set(os_dict[next_rxn]['reactants'].keys())
                    else:
                        next_reacts = set(os_dict[next_rxn]['products'].keys())
                    difference = set(cur_prods.keys()).intersection(set(cur_prods.keys()).symmetric_difference(next_reacts))
                else:
                    difference = set(cur_prods.keys()).intersection(set(cur_prods.keys()).symmetric_difference(set([self.target])))

                for cpd_id in difference:
                    cpd = cur_prods[cpd_id]
                    if cpd == 'None':
                        cpd_name = re.sub('_\w+\d+$', '', cpd)
                        cpd_name = '_'.join(cpd_name.split('/'))
                        if cpd_name not in promiscuous_products[rxn]:
                            misc_products[rxn] += [cpd_id]
                    elif cpd not in promiscuous_products[rxn]:
                        misc_products[rxn] += [cpd]

            get_info(path.root)

            if _FBA:
                max_flux = max(fba_values.values())
                if max_flux != 0:
                    fba_values = {fv_key: float(fv_val/max_flux) for fv_key, fv_val in fba_values.iteritems()}
            
            
            last_id=0
            cdxml_compounds_path = self.output_path + '/compounds/'
            def generate_cdxml(root, last_id=last_id):
                previous = []
                for child in root.children:
                    child_box, last_id = generate_cdxml(child, last_id=last_id)
                    previous += [child_box]

                rxn = root.root
                CE = CDXML_Editor(cdxml_files_path=cdxml_compounds_path)
                last_id = CE.add_reactants(cdxml_reactants[rxn], previous, last_id)

                if reaction_proteins[rxn]:
                    CE.add_transition(promiscuous_reactants[rxn],promiscuous_products[rxn],misc_products[rxn], 
                                    reaction_proteins=reaction_proteins[rxn])
                elif show_rxn_info:
                    CE.add_transition(promiscuous_reactants[rxn],promiscuous_products[rxn],misc_products[rxn], 
                                    reaction_solvents=reaction_solvents[rxn],reaction_catalysts=reaction_catalysts[rxn],reaction_SPRESI_info=reaction_SPRESI_info[rxn])
                else:
                    CE.add_transition(promiscuous_reactants[rxn],promiscuous_products[rxn],misc_products[rxn])
                
                if _FBA:
                    MAIN_CE.add_color(fba_values[rxn])
                    CE.set_FBA(MAIN_CE.color_index)
                
                CE.set_products(cdxml_products[rxn])                
                return CE, last_id

            pathway_cdxml, last_id = generate_cdxml(path.root)

            # FINALLY ADD TARGET
            target_CE = CDXML_Editor(cdxml_files_path=cdxml_compounds_path)
            target_name = self.DB.get_compound_name(self.target)
            if target_name == 'None':
                target_name = '_'.join(self.target.split('/'))            
            target_CE.add_product(target_name, last_id)
            pathway_cdxml.append(target_CE.container, arrange="right")

            MAIN_CE.append(pathway_cdxml.container)
            MAIN_CE.generate_file()
        

    # def generate_reaction_SMILES(self):
    #     '''Produce files that have reaction smiles for each reaction in a pathway'''
    #     ordered_paths = self.order_of_paths()
    #     self.ordered_paths = ordered_paths
    #     if self.figures:
    #         target_reformat_orig = re.sub('/', '_', self.target)
    #         self.generate_output_folders(self.output_path+'/solution_smiles')
    #         self.generate_output_folders(self.output_path+'/solution_cdxml')
    #         self.generate_output_folders(self.output_path+'/solution_cdxml/'+
    #                                      target_reformat_orig+'_solutions')
    #         self.generate_output_folders(self.output_path+'/solution_figures')
    #         self.generate_output_folders(self.output_path+'/solution_figures/'+
    #                                      target_reformat_orig+'_solutions')
            
    #         target_reformat = self.alter_name_length(self.output_path+'/solution_smiles/reaction_smile_'+
    #                                                  target_reformat_orig+'_solution_1.smi',
    #                                                  target_reformat_orig)
    #         for count_pathway, os_dict in self.reactions.iteritems():
    #             smiles_reactants = {}
    #             smiles_products = {}

    #             with open(self.output_path+'/solution_smiles/reaction_smile_'+
    #                       target_reformat+'_solution_'+str(count_pathway)+
    #                       '.smi', 'w') as fout:
    #                 array = self.IN.createArray()
    #                 for counter in reversed(ordered_paths[count_pathway].keys()):
    #                     rxn = ordered_paths[count_pathway][counter]
    #                     smiles_reactants[rxn] = []
    #                     smiles_products[rxn] = []
    #                     if os_dict[rxn]['direction'] == 'forward':
    #                         smiles_reactants[rxn] = self.get_SMILES(os_dict[rxn]['reactants'],
    #                                                                 smiles_reactants[rxn])
    #                         smiles_products[rxn] = self.get_SMILES(os_dict[rxn]['products'],
    #                                                                smiles_products[rxn])
    #                     else:
    #                         smiles_reactants[rxn] = self.get_SMILES(os_dict[rxn]['products'],
    #                                                                 smiles_reactants[rxn])
    #                         smiles_products[rxn] = self.get_SMILES(os_dict[rxn]['reactants'],
    #                                                                smiles_products[rxn])
    #                 for count_rxn in reversed(ordered_paths[count_pathway].keys()):
    #                     key = ordered_paths[count_pathway][count_rxn]
    #                     fout.write('.'.join(smiles_reactants[key])+'>>'+
    #                                '.'.join(smiles_products[key])+'\n')
    #                     rxn = self.IN.loadReaction('.'.join(smiles_reactants[key])+'>>'+
    #                                                '.'.join(smiles_products[key]))
    #                     self.IN.setOption('render-grid-title-property', str(key)+
    #                                       ' '+str(self.DB.get_reaction_name(key)))
    #                     array.arrayAdd(rxn)
    #                     self.IR.renderToFile(rxn, self.output_path+'/solution_cdxml/'+
    #                                          target_reformat_orig+'_solutions/'+
    #                                          'Solution_'+str(count_pathway)+'_rxn_'+
    #                                          str(count_rxn)+'_'+self.target_organism_name+'.cdxml')

    #                 self.IN.setOption("render-output-format", "cdxml")
    #                 self.IR.renderGridToFile(array, None, len(smiles_reactants.keys()),
    #                                          self.output_path+'/solution_figures/'+target_reformat_orig+
    #                                          '_solutions/'+'Solution_'+str(count_pathway)+'_'+
    #                                          self.target_organism_name+'.cdxml')