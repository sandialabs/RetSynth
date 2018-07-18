from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Check based on compound structure if it is in the database'

import re
from tqdm import tqdm
from copy import deepcopy
from sys import platform
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_inchi

def verbose_print(verbose, line):
    if verbose:
        print(line)

class TanimotoStructureSimilarity(object):
    """Identifies if compounds are in the database but under a different ID
    based on structure similarity"""
    def __init__(self, targets, all_compounds, cytosol, extracellular, verbose, threshold_score = 1):
        """Initialize"""
        self.verbose = verbose
        self.targets = []
        self.track_final_cpd = []
        '''Remove duplicated targets'''
        for target in targets:
            if target not in self.targets:
                self.targets.append(target)
        self.individualtargets = []
        self.organisms = []
        for target in self.targets:
            self.individualtargets.append(target[0])
            organism = ','.join([target[1], target[2], target[3]])
            if organism not in self.organisms:
                self.organisms.append(organism)
        self.organisms = list(set(self.organisms))
        self.all_compounds = all_compounds
        self.cytosol = cytosol
        self.extracellular = extracellular
        self.threshold_score = threshold_score
        self.IN = indigo.Indigo()
        self.INCHI = indigo_inchi.IndigoInchi(self.IN)
        self.calculate_tanimoto_score()

    def remove_compartment_info_from_cpdID(self, cpd):
        cpd = re.sub('\_\w{1}0$', '', cpd)
        cpd = re.sub('CCO__45__MIT__45__IM__45__SPC__45__CCO__45__MIT__45__LUM', '', cpd)
        cpd = re.sub('CCO__45__CYTOSOL__45__CCO__45__CHLOR__45__STR', '', cpd)
        cpd = re.sub('_CCO__45__CW__45__BAC__45__POS0', '', cpd)
        cpd = re.sub('CCI__45__PERI__45__BAC__45__GN', '', cpd)
        cpd = re.sub('CCO__45__OUT__45__CCO__45__IN', '', cpd)
        cpd = re.sub('CCO__45__PM__45__BAC__45__NEG', '', cpd)
        cpd = re.sub('CCO__45__CW__45__BAC__45__POS', '', cpd)
        cpd = re.sub('CCO__45__RGH__45__ER__45__LUM', '', cpd)
        cpd = re.sub('_CCO__45__PLASTID__45__STR0', '', cpd)
        cpd = re.sub('CCO__45__UNKNOWN__45__SPACE', '', cpd)
        cpd = re.sub('_CCO__45__PEROX__45__MEM0', '', cpd)
        cpd = re.sub('CCO__45__PLASTID__45__STR', '', cpd)
        cpd = re.sub('CCO__45__CHROM__45__STR', '', cpd)
        cpd = re.sub('CCO__45__CHROM__45__STR', '', cpd)
        cpd = re.sub('CCO__45__GOLGI__45__LUM', '', cpd)
        cpd = re.sub('CCO__45__PEROX__45__MEM', '', cpd)
        cpd = re.sub('CCO__45__MIT__45__IMEM', '', cpd)
        cpd = re.sub('CCO__45__VESICLE', '', cpd)
        cpd = re.sub('CCO-PEROX-LUM', '', cpd)
        cpd = re.sub('_CCO__45__IN0', '', cpd)
        cpd = re.sub('CCO__45__OUT', '', cpd)
        cpd = re.sub('CCO__45__IN', '', cpd)
        cpd = re.sub('\_0$', '', cpd)
        return cpd

    def process_targets(self, targets):
        """Reformat targets to remove compartment information"""
        reformat_target = []
        for target in targets:
            target = self.remove_compartment_info_from_cpdID(target)
            reformat_target.append(target)
        return(reformat_target)

    def get_original_target(self, tmol):
        """Obtain original target"""
        index = None
        for count, target in enumerate(self.targets):
            if tmol in target:
                index =  count
         
        if not index: 
            verbose_print(self.verbose, 'Could not get target '+tmol)
            return (None)
        else:
            return index

    def get_tanimoto_score(self, tmol, threshold):
        temp = {}
        for db_cpd in self.db_cpds_fp:
            score = self.IN.similarity(self.individualtargets_p_fp[tmol], self.db_cpds_fp[db_cpd], 'tanimoto')
            temp[db_cpd] = score
        max_score_cpds = [temp.keys()[temp.values().index(i)] for i in temp.values() if float(i) >= float(threshold)]
        return (max_score_cpds)

    def calculate_tanimoto_score(self):
        """Calculate similarity for compounds that are not in the database"""
           
        self.finaltargets = self.targets
        self.individualtargets_p = self.process_targets(set(self.individualtargets))
        self.retrieve_fingerprints()
        for count ,tmol in enumerate(self.individualtargets_p_fp.keys()):
            if tmol+'_'+self.cytosol in self.all_compounds:
                max_score_cpds = self.get_tanimoto_score(tmol, 1)
                if max_score_cpds:
                    for max_score_cpd in set(max_score_cpds):
                        if max_score_cpd != tmol:
                            new_target = self.extract_db_cpd_ID(max_score_cpd)
                            if new_target not in self.individualtargets and new_target not in self.track_final_cpd:
                                verbose_print(self.verbose, "STATUS: The following compound {} have 100 percent similarity with the original target {} therefore adding them to to targets".format(max_score_cpd, tmol))                
                                self.fill_final_targets(new_target)
            
            elif tmol in self.cpd_dict.keys():
                index = self.get_original_target(tmol+'_'+self.cytosol)
                if index:
                    new_target = self.extract_db_cpd_ID(tmol)
                    verbose_print(self.verbose, "STATUS: Switching target ID from {} to {}".format(tmol+'_'+self.cytosol, new_target))
                    self.finaltargets[index][0] = new_target

            else:
                max_score_cpds = self.get_tanimoto_score(tmol, self.threshold_score)
                if max_score_cpds:
                    index = self.get_original_target(tmol+'_'+self.cytosol)
                    if index:
                        del self.finaltargets[index]
                    verbose_print(self.verbose, 'STATUS: {} compounds have {} or greater similarity to target compound {}'.format(len(set(max_score_cpds)), float(self.threshold_score)*100, tmol+'_'+self.cytosol))
                    for max_score_cpd in set(max_score_cpds):
                        new_target = self.extract_db_cpd_ID(max_score_cpd)
                        if new_target not in self.individualtargets and new_target not in self.track_final_cpd:
                            print ('STATUS: Adding compound {} to target list'.format(max_score_cpd))
                            self.fill_final_targets(new_target)
                else:
                    verbose_print(self.verbose,'STATUS: No compounds in the database that are {} percent similar to target {}'.format(float(self.threshold_score)*100, tmol+'_'+self.cytosol))

    def fill_final_targets(self, new_target):
        for organism in self.organisms:
            organisms = organism.split(',')
            self.finaltargets.append([new_target]+organisms)
            self.track_final_cpd.append(new_target)

    def retrieve_fingerprints(self):
        """Retrieve fingerprints for database compounds and targets"""
        db_cpds = []
        self.cpd_dict = {}
        for cpd in self.all_compounds:
            originalcpd = deepcopy(cpd)
            cpd = self.remove_compartment_info_from_cpdID(cpd)
            db_cpds.append(cpd)
            self.cpd_dict.setdefault(cpd, []).append(originalcpd)

        self.individualtargets_p_fp = {}
        print ("STATUS: getting fingerprints for target compounds")
        for count, tmol in enumerate(tqdm(self.individualtargets_p)):
            try:
                self.individualtargets_p_fp[tmol] = self.INCHI.loadMolecule(tmol).fingerprint('full')
            except indigo.IndigoException:
                verbose_print(self.verbose, 'Could not get fingerprint for {}'.format(tmol))

        self.db_cpds_fp = {}
        db_cpds_set = set(db_cpds)
        print ("STATUS: getting fingerprints for database compounds")        
        for db_cpd in tqdm(db_cpds_set):
            try:
                mol = self.INCHI.loadMolecule(db_cpd)
                self.db_cpds_fp[db_cpd] = mol.fingerprint('full')
            except indigo.IndigoException:
                pass

    def extract_db_cpd_ID(self, cpdID):
        """Extract correct database ID for a compound"""
        if cpdID+'_'+self.cytosol in self.cpd_dict[cpdID]:
            new_target = cpdID+'_'+self.cytosol
        else:
            if len(self.cpd_dict[cpdID]) > 1:
                for value in self.cpd_dict[cpdID]:
                    if value.endswith(self.extracellular):
                        pass
                    else:
                        new_target = value
            else:
                new_target = self.cpd_dict[cpdID][0]
        return(new_target)
