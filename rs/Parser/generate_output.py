from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Generate output'

import sys
import re
import os
import csv
import openpyxl
import shutil
csv.field_size_limit(sys.maxsize)

class Output(object):
    """Opens and fills output files produced by software"""
    def __init__(self, db, output_path, FBA=False, KO=False, raw_solutions=True):
        '''Initialize class: generates new output files for this analysis of rs'''
        self.DB = db
        self.FBA = FBA
        self.KO = KO
        self.output_path = output_path
        self.optimal_paths = open(output_path+'/optimal_pathways.txt', 'w')
        self.optimal_paths.close()
        if self.FBA:
            self.activemetabolism = open(output_path+'/active_metabolism.txt', 'w')
            self.flux_ouptput = open(output_path+'/flux_output.txt', 'w')
            self.flux_individual_output = open(output_path+'/flux_individualfluxes_output.txt', 'w')
            self.theoyield = open(output_path+'/theoretical_yield.txt', 'w')

            self.activemetabolism.close()
            self.flux_ouptput.close()
            self.flux_individual_output.close()
            self.theoyield.close()

        if self.KO:
            self.essentialrxns = open(output_path+'/essentialrxns_output.txt', 'w')
            self.fluxKO_ouptput = open(output_path+'/fluxKO_output.txt', 'w')
            self.fluxKO_ty_output = open(output_path+'/fluxKO_theoreticalyields_output.txt', 'w')
            self.fluxKO_in_ty_output = open(output_path+'/fluxKO_increased_theoreticalyields_output.txt', 'w')
            self.fluxKO_in_ty_output.close()
            self.essentialrxns.close()
            self.fluxKO_ty_output.close()
            self.fluxKO_ouptput.close()

        if raw_solutions:
            try:
                os.mkdir(output_path+'/raw_compound_solutions')
            except OSError:
                shutil.rmtree(output_path+'/raw_compound_solutions')
                os.mkdir(output_path+'/raw_compound_solutions')

    def output_open_paths_all_organism_file(self, target_compound_ID):
        '''
        Open file to generate file that summarizes how many reactions would
        need to be added to every organism in the daabase to produce a target organism
        '''
        file_name = 'path_length_all_organism_'+target_compound_ID+'.txt'
        self.all_organisms = open(self.output_path+'/'+file_name, 'w')
        self.all_organisms.close()

    def output_raw_solutions(self, compound, org_ID, ordered_paths, reactions, incpds):
        cpdname = self.DB.get_compound_name(compound)
        org_name = self.DB.get_organism_name(org_ID)
        if cpdname == 'None' or cpdname == 'none' or cpdname == '':
            cpdname = re.sub('/', '_', compound)

        def add_info_to_output(reactants, products, fin, protein, incpds):
            '''adding reaction info to output file '''
            for react in reactants:
                type_cpd = 'Ext'
                if react in incpds:
                    type_cpd = 'Int'
                if protein != 'None':
                    line = ', '.join([react, protein, str(len(os_dict)), str(count_pathway), type_cpd])
                else:
                    line = ', '.join([react, rxn, str(len(os_dict)), str(count_pathway), type_cpd])
                fin.write(line+'\n')

            for prod in products:
                type_cpd = 'Ext'
                if prod in incpds:
                    type_cpd = 'Int'
                if protein != 'None':
                    line = ', '.join([protein, prod, str(len(os_dict)), str(count_pathway), type_cpd])
                else:
                    line = ', '.join([rxn, prod, str(len(os_dict)), str(count_pathway),  type_cpd])
                fin.write(line+'\n')
        
        def alter_name_length(path_to_figure, cpdname):
            '''Shorten compound name if it is too long'''
            if len(path_to_figure) > 250:
                remove_variable = len(path_to_figure) - 250
                cpdname = cpdname[:-remove_variable]
            return cpdname

        cpdname = alter_name_length(self.output_path+'/raw_compound_solutions/compound_'+cpdname+'_'+org_name+'_outputfile.txt', cpdname)
        with open(self.output_path+'/raw_compound_solutions/compound_'+cpdname+'_'+org_name+'_outputfile.txt', 'w') as fin:
            for count_pathway, os_dict in reactions.iteritems():
                for counter in reversed(ordered_paths[count_pathway].keys()):
                    rxn = ordered_paths[count_pathway][counter]
                    org = os_dict[rxn]['organisms'][0]
                    protein = self.DB.get_proteins(rxn, org)
                    if os_dict[rxn]['direction'] == 'forward':
                        add_info_to_output(os_dict[rxn]['reactants'], os_dict[rxn]['products'], fin, protein, incpds)

                    else:
                        add_info_to_output(os_dict[rxn]['products'], os_dict[rxn]['reactants'], fin, protein, incpds)

    def output_paths_all_organisms(self, target_compound_ID, pathlength, orgs, org_names):
        '''
        Outputs information on the number of reactions that need to be added to an
        organism to get a target compound, this file is only generated all organisms
        are being examined to see if they can produce target compound
        '''
        with open(self.output_path+'/'+'path_length_all_organism_'+target_compound_ID+'.txt', 'a') as self.all_organisms:
            self.all_organisms.write('{} reaction steps need to be added to get {} in organism {} ({})\n'.format(pathlength,
                                                                                                                 target_compound_ID,
                                                                                                                 ','.join(orgs),
                                                                                                                 ','.join(org_names)))
    def output_compound_natively_present_in_target_organism(self, target_info):
        '''
        Outputs information if a target compound is already present in an organism
        '''
        with open(self.output_path+'/optimal_pathways.txt', 'a') as self.optimal_paths:
            print ('{} in species {} already'.format(target_info[0], target_info[2]))
            self.optimal_paths.write('{} in species {} already\n'.format(target_info[0], target_info[2]))

    def output_shortest_paths(self, target_info, temp_rxns, tox_excpds=None):
        #self.optimal_paths = optimal_paths
        # self.optimal_paths = open(self.output_path+'/optimal_pathways.txt', 'a')
        '''
        Outputs reactions and compounds that need to be
        added to an organism to get target compound
        '''
        t = target_info[0]
        target_org = target_info[2]
        with open(self.output_path+'/optimal_pathways.txt', 'a') as self.optimal_paths:
            if len(temp_rxns) == 0:
                print ('No paths could be found to get to target compound {} {} in target organism {}'.format(t,
                                                                                                              self.DB.get_compound_name(t),
                                                                                                              self.DB.get_organism_name(target_org)))

                self.optimal_paths.write('No paths could be found to get to target compound {} {} in target organism {}\n'.format(t,
                                                                                                                                  self.DB.get_compound_name(t),
                                                                                                                                  self.DB.get_organism_name(target_org)))
            else:
                print ('\nSHORTEST PATH FOR {} {} in target organism {}'.format(t, self.DB.get_compound_name(t),
                                                                                self.DB.get_organism_name(target_org)))
                self.optimal_paths.write('\nSHORTEST PATH FOR {} {} in target organism {}\n'.format(t, self.DB.get_compound_name(t),
                                                                                                    self.DB.get_organism_name(target_org)))
                for count, os_dict in temp_rxns.iteritems():
                    print ('Solution {}'.format(count))
                    self.optimal_paths.write('Solution {}\n'.format(count))
                    for r in os_dict:
                        if r.endswith('_s'):
                            print ('\t'.join([r, os_dict[r]['name'], os_dict[r]['direction']]))
                            self.optimal_paths.write('\t'.join([r, os_dict[r]['name'],
                                                                os_dict[r]['direction'],
                                                                ','.join(self.DB.get_solvents(r))+':solvents',
                                                                ','.join(self.DB.get_catalysts(r))+':catalysts',
                                                                ','.join(self.DB.get_time(r))+':time',
                                                                ','.join(self.DB.get_temperature(r))+':temperature',
                                                                ','.join(self.DB.get_pressure(r))+':pressure',
                                                                ','.join(self.DB.get_yield(r))+':yield',
                                                                ','.join(self.DB.get_reference(r))+':reference'])+ '\n')
                        else:
                            proteins = self.DB.get_proteins(r, os_dict[r]['organisms'][0])
                            proteins = re.sub('\(', '', proteins)
                            proteins = re.sub('\)', '', proteins)
                            proteinslist = proteins.split(' ')
                            finalproteinlist = []
                            for protein in proteinslist:
                                finalproteinlist.append(protein)
                            finalprotein = ' '.join(finalproteinlist)

                            genes = self.DB.get_genes(r, os_dict[r]['organisms'][0])
                            genes = re.sub('\(', '', genes)
                            genes = re.sub('\)', '', genes)
                            geneslist = genes.split(' ')
                            finalgenelist = []
                            for gene in geneslist:
                                finalgenelist.append(gene)
                            finalgene = ' '.join(finalgenelist)
                            print ('\t'.join([r, os_dict[r]['name'], os_dict[r]['direction'], finalprotein, finalgene]))
                            self.optimal_paths.write('\t'.join([r, os_dict[r]['name'],
                                                                os_dict[r]['direction'], finalprotein, finalgene,
                                                                str(len(os_dict[r]['organisms']))+
                                                                ' number of species that contain this reaction',
                                                                ','.join(os_dict[r]['organisms'])])+'\n')

                        if tox_excpds:
                            for react in os_dict[r]['reactants']:
                                try:
                                    print ('\t{}\t{} reactants\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}'.format(react, os_dict[r]['reactants'][react],
                                                                                                                                                   tox_excpds[react][0], tox_excpds[react][1][0], tox_excpds[react][1][1]))
                                    self.optimal_paths.write('\t{}\t{} reactants\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}\n'.format(react,
                                                                                                                                                                       os_dict[r]['reactants'][react]),
                                                                                                                                                                       tox_excpds[react][0],
                                                                                                                                                                       tox_excpds[react][1][0],
                                                                                                                                                                       tox_excpds[react][1][1])
                                except IndexError:
                                    print ('\t{}\t{} reactants\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}'.format(react, os_dict[r]['reactants'][react],
                                                                                                                                                   'None', 'None', 'None'))
                                    self.optimal_paths.write('\t{}\t{} reactants\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}\n'.format(react,
                                                                                                                                                                       os_dict[r]['reactants'][react],
                                                                                                                                                                       'None', 'None', 'None'))
                            for prod in os_dict[r]['products']:
                                try:
                                    print ('\t{}\t{} products\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}'.format(prod, os_dict[r]['products'][prod],
                                                                                                                                                  tox_excpds[prod][0], tox_excpds[prod][1][0], tox_excpds[prod][1][1]))
                                    self.optimal_paths.write('\t{}\t{} products\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}\n'.format(prod, os_dict[r]['products'][prod],
                                                                                                                                                                      tox_excpds[prod][0], tox_excpds[prod][1][0],
                                                                                                                                                                      tox_excpds[prod][1][1]))
                                except IndexError:
                                    print ('\t{}\t{} products\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}'.format(prod, os_dict[r]['products'][prod],
                                                                                                                                                  'None', 'None', 'None'))
                                    self.optimal_paths.write('\t{}\t{} products\ttoxicity = classification-{}, probability-High toxicity {} Low toxicity {}\n'.format(prod, os_dict[r]['products'][prod],
                                                                                                                                                                      'None', 'None', 'None'))
                        else:
                            for react in os_dict[r]['reactants']:
                                print ('\t{}\t{} reactants'.format(react, os_dict[r]['reactants'][react]))
                                self.optimal_paths.write('\t{}\t{} reactants\n'.format(react,
                                                                                       os_dict[r]['reactants'][react]))
                            for prod in os_dict[r]['products']:
                                print ('\t{}\t{} products'.format(prod, os_dict[r]['products'][prod]))
                                self.optimal_paths.write('\t{}\t{} products\n'.format(prod, os_dict[r]['products'][prod]))

    def output_activemetabolism(self, am):
        '''
        Outputs metabolites that are actively produced in a specific compund
        '''
        with open(self.output_path+'/active_metabolism.txt', 'a') as self.activemetabolism:
            for key, values in am.iteritems():
                self.activemetabolism.write('Active compounds in {}\n'.format(key))
                for cpd in values[0]:
                    self.activemetabolism.write('{}\t{}\n'.format(cpd, self.DB.get_compound_name(cpd)))


    def output_FBA(self, target_info, org_fbasolution, optimized_fba, comparisonresults, temp):
        '''
        Outputs reactions and corresponding fluxes for fluxes that change significantly
        between wild-type FBA simulation (simulation without added reactions) and
        FBA simulation with added reactions to produce target compound
        '''
        if target_info[1] != '':
            target = target_info[1]
        else:
            target = target_info[0]
        with open(self.output_path+'/flux_individualfluxes_output.txt', 'a') as self.flux_individual_output:
            print ('FBA Solutions for {}'.format(target))
            print ('{}\t{} objective function solutions for wild-type and mutant'.format(round(org_fbasolution.f, 2),
                                                                                         round(optimized_fba.fbasol.f, 2)))
            self.flux_individual_output.write('FBA Solutions for {}'.format(target)+'\n')
            self.flux_individual_output.write('{}\t{} objective function solutions for wild-type and mutant\n'.format(round(org_fbasolution.f, 2),
                                                                                                                      round(optimized_fba.fbasol.f, 2)))
            for x, value in optimized_fba.fbasol.x_dict.iteritems():
                self.flux_individual_output.write('{}\t{}\n'.format(x, value))

        with open(self.output_path+'/flux_output.txt', 'a') as self.flux_ouptput:
            self.flux_ouptput.write('FBA Solutions for {}\n'.format(target))
            self.flux_ouptput.write('{}\t{} objective function solutions for wild-type and mutant\n'.format(round(org_fbasolution.f, 2),
                                                                                                            round(optimized_fba.fbasol.f, 2)))
            self.flux_ouptput.write('\nFluxes that differ by 1.5 fold for reactions between wildtype and mutant:\n')
            self.flux_ouptput.write('\t\t\twildtype flux\tmutantflux\n')
            for x, fluxvalue in comparisonresults.fluxchange.iteritems():
                self.flux_ouptput.write('\t{}\t{}\n'.format(x, fluxvalue))
            self.flux_ouptput.write('\nFluxes for added reactions in mutant:'+'\n')
            for r, value in comparisonresults.externalrxnfluxes.iteritems():
                self.flux_ouptput.write('\t{}\t{}\n'.format(r, value))
            if comparisonresults.maxpath == 'No added path':
                print ('{} - compound could be produced in target organism'.format(comparisonresults.maxpath))
                print ('\t {} production of target compound'.format(comparisonresults.maxflux))
                self.flux_ouptput.write('\t {} - compound could be produced in target organism\n'.format(comparisonresults.maxpath))
                self.flux_ouptput.write('\t {} production of target compound\n'.format(comparisonresults.maxflux))
            else:
                self.flux_ouptput.write('\nExternal pathway with most flux:\n')
                self.flux_ouptput.write('\tPath {}\t{}\n\tTotal flux through path: {}\n'.format(comparisonresults.maxpath,
                                                                                                temp[comparisonresults.maxpath].keys(),
                                                                                                comparisonresults.maxflux))
                print ('External pathway with most flux:')
                print ('\tPath {}\t{}\n\tTotal flux through path: {}'.format(comparisonresults.maxpath,
                                                                             temp[comparisonresults.maxpath].keys(),
                                                                             comparisonresults.maxflux))

    def output_FBA_KOs(self, target_info, fbasolution, compound_dict, comparisonKOresults, temp):
        '''
        When reaction knockouts are performed list all reaction that have a significant
        flux difference from the orignal FBA simulation
        '''
        if target_info[1] != '':
            target = target_info[1]
        else:
            target = target_info[0]
        with open(self.output_path+'/fluxKO_output.txt', 'a') as self.fluxKO_ouptput:
            self.fluxKO_ouptput.write('{} target compound\n'.format(target))
            self.fluxKO_ouptput.write('Fluxes that differ by 1.5 fold for reactions between wildtype and mutant:\n')
            for r, value in comparisonKOresults.fluxchange.iteritems():
                self.fluxKO_ouptput.write('\t {} knockout\n'.format(r))
                self.fluxKO_ouptput.write('\t\t\t\twildtype flux\tmutantflux\n')
                for rk in value:
                    self.fluxKO_ouptput.write('\t\t{}\t{}\n'.format(rk, value[rk]))
                if comparisonKOresults.maxpath[r] == 'No added path':
                    self.fluxKO_ouptput.write('\t\t{} - compound could be produced in target organism\n'.format(comparisonKOresults.maxpath[r]))
                    self.fluxKO_ouptput.write('\t\t{} production of target compound\n'.format(comparisonKOresults.maxflux[r]))
                else:
                    self.fluxKO_ouptput.write('\t\tPath {}\t{}\n'.format(comparisonKOresults.maxpath[r],
                                                                         temp[comparisonKOresults.maxpath[r]].keys()))
                    self.fluxKO_ouptput.write('\t\t'+'Total flux through path: {}\n'.format(comparisonKOresults.maxflux[r]))

        objectivesol = fbasolution.x_dict['Sink_'+compound_dict[target]]
        glucose = True
        try:
            glucoseimport = fbasolution.x_dict['EX_cpd00027_e0']
        except KeyError:
            print ('No transport reaction glucose rxn')
            glucose = False
        if glucose:
            try:
                wt_ty = abs(round(round(objectivesol, 2)/round(glucoseimport,2), 2))
            except ZeroDivisionError:
                wt_ty = 'NA'
        else:
            wt_ty = 'NA'

        with open(self.output_path+'/fluxKO_theoreticalyields_output.txt') as self.fluxKO_ty_output:
            line = self.fluxKO_ty_output.readline()

        with open(self.output_path+'/fluxKO_theoreticalyields_output.txt', 'a') as self.fluxKO_ty_output:
            if line.startswith('#'):
                self.ko_ty = '\t'.join(format(x, "10.3f") for x in comparisonKOresults.objective_function_ko.values())
                self.fluxKO_ty_output.write(target+'-'+self.DB.get_compound_name(target)+'\t'+target_info[2]+'-'+self.DB.get_organism_name(target_info[2])+'\t'+str(wt_ty)+'\t'+self.ko_ty+'\n')
            else:
                self.fluxKO_ty_output.write('#Target ID\tOrganism ID\twild type theoretical yield\t'+'\t'.join(comparisonKOresults.objective_function_ko.keys())+'\n')
                self.ko_ty = '\t'.join(format(x, "10.3f") for x in comparisonKOresults.objective_function_ko.values())
                self.fluxKO_ty_output.write(target+'-'+self.DB.get_compound_name(target)+'\t'+target_info[2]+'-'+self.DB.get_organism_name(target_info[2])+'\t'+str(wt_ty)+'\t'+self.ko_ty+'\n')                

        with open(self.output_path+'/fluxKO_increased_theoreticalyields_output.txt', 'a') as self.fluxKO_in_ty_output:
            if wt_ty != 'NA':
                count=0
                for rko, value in comparisonKOresults.objective_function_ko.iteritems():
                    if rko.startswith('EX'):
                        pass
                    else:
                        if value > wt_ty:
                            count+=1
                            if count == 1:
                                self.fluxKO_in_ty_output.write('Reaction knockouts that increase theoretical yield of {} of compound {} in organism {}\n'.format(wt_ty, 
                                                                                                                                                             target+'-'+self.DB.get_compound_name(target),
                                                                                                                                                             target_info[2]+'-'+self.DB.get_organism_name(target_info[2])))
                                self.fluxKO_in_ty_output.write("reaction knockout ID\treaction knockout name\tcatalytic genes\tcatalytic proteins\treaction formula\tyield\n")
                            self.fluxKO_in_ty_output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rko, self.DB.get_reaction_name(rko),
                                                                                             self.DB.get_genes(rko, target_info[2]),
                                                                                             self.DB.get_proteins(rko, target_info[2]),
                                                                                             ','.join(self.DB.get_reactants(rko))+'-->'+','.join(self.DB.get_products(rko)),
                                                                                             value))
 
    def output_essential_reactions(self, target_compound_ID, target_organism_ID, er):
        '''
        When reaction knockouts are performed, outputs all reactions that when removed
        cause decrease in target production
        '''
        with open(self.output_path+'/essentialrxns_output.txt', 'a') as self.essentialrxns:
            self.essentialrxns.write('Essential rxns for production of {} in {}\n'.format(target_compound_ID,
                                                                                          target_organism_ID))
            for rxn in er:
                if rxn.startswith('EX'):
                    pass
                else:
                    self.essentialrxns.write('{}\t{}\n'.format(rxn, self.DB.get_reaction_name(rxn)))

    def output_theoretical_yield(self, target_compound_ID, target_organism_ID,
                                 fbasolution, compounds_dict):
        '''
        Generates output file with theoretical yields
        '''
        objectivesol = fbasolution.x_dict['Sink_'+compounds_dict[target_compound_ID]]
        glucose = True
        ##########CHECK FOR GLUCOSE IMPORT REACTION (CURRENTLY KBASE RXN))###########
        try:
            glucoseimport = fbasolution.x_dict['EX_cpd00027_e0']
        except KeyError:
            print ('No transport reaction glucose rxn')
            glucose = False
        ##########CHECK FOR BIOMASS REACTION (CURRENTLY KBASE RXN))###########
        try:
            biomassrxn = fbasolution.x_dict['biomass0_'+target_organism_ID]
        except KeyError:
            biomassrxn = 'NA'
        ##########CALCULATE WT THEORETICAL YIELD###########
        if glucose:
            try:
                wt_ty = abs(round(round(objectivesol, 2)/round(glucoseimport,2), 2))
            except ZeroDivisionError:
                wt_ty = 'NA'
        else: 
            wt_ty = 'NA'

        with open(self.output_path+'/theoretical_yield.txt', 'a') as self.theoyield:
            self.theoyield.write('{}-{}\t{}-{}\tGlucose Flux: {}\tTarget Production: {}\tTheoretical Yield: {} mol {} /mol glucose\tBiomass: {}\n'.format(target_compound_ID,
                                                                                                                                                          self.DB.get_compound_name(target_compound_ID),                                                                                                                                                         target_organism_ID,
                                                                                                                                                          self.DB.get_organism_name(target_organism_ID),
                                                                                                                                                          round(glucoseimport, 2),
                                                                                                                                                          round(objectivesol, 2),
                                                                                                                                                          wt_ty, target_compound_ID, biomassrxn))
    def convert_output_2_xlsx(self):
        '''converts txt files to output files'''
        print ('STATUS: Converting output text files to xlsx format')
        def convert_output_txt_files_2_xlsx(input_file, output_file):
            wb = openpyxl.Workbook()
            ws = wb.worksheets[0]

            with open(input_file, 'rb') as data:
                reader = csv.reader(data, delimiter='\t')
                for row in reader:
                    ws.append(row)
            wb.save(output_file)

        convert_output_txt_files_2_xlsx(self.output_path+'/optimal_pathways.txt', self.output_path+'/optimal_pathways.xlsx')
        if self.FBA:
            convert_output_txt_files_2_xlsx(self.output_path+'/active_metabolism.txt', self.output_path+'/active_metabolism.xlsx')
            convert_output_txt_files_2_xlsx(self.output_path+'/flux_individualfluxes_output.txt', self.output_path+'/flux_individualfluxes_output.xlsx')
            convert_output_txt_files_2_xlsx(self.output_path+'/flux_output.txt', self.output_path+'/flux_output.xlsx')
            convert_output_txt_files_2_xlsx(self.output_path+'/theoretical_yield.txt', self.output_path+'/theoretical_yield.xlsx')
        if self.KO:
            convert_output_txt_files_2_xlsx(self.output_path+'/essentialrxns_output.txt', self.output_path+'/essentialrxns_output.xlsx')
            convert_output_txt_files_2_xlsx(self.output_path+'/fluxKO_output.txt', self.output_path+'/fluxKO_output.xlsx')
            convert_output_txt_files_2_xlsx(self.output_path+'/fluxKO_theoreticalyields_output.txt', self.output_path+'/fluxKO_theoreticalyields_output.xlsx')
            convert_output_txt_files_2_xlsx(self.output_path+'/fluxKO_increased_theoreticalyields_output.txt', self.output_path+'/fluxKO_increased_theoreticalyields_output.xlsx')

        