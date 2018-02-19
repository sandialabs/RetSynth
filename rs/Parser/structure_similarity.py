from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Check based on compound structure if it is in the database'

import re
from copy import deepcopy
from sys import platform
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_inchi


class TanimotoStructureSimilarity(object):
    """Identifies if compounds are in the database but under a different ID
    based on structure similarity"""
    def __init__(self, targets, all_compounds, cytosol, extracellular, threshold_score = 1):
        """Initialize"""
        self.targets = []
        '''Remove duplicated targets'''
        for target in targets:
            if target not in self.targets:
                self.targets.append(target)
        self.individualtargets = []
        for target in self.targets:
            self.individualtargets.append(target[0])
        self.all_compounds = all_compounds
        self.cytosol = cytosol
        self.extracellular = extracellular
        self.threshold_score = threshold_score
        self.IN = indigo.Indigo()
        self.INCHI = indigo_inchi.IndigoInchi(self.IN)
        self.calculate_tanimoto_score()

    def process_targets(self, targets):
        """Reformat targets to remove compartment information"""
        reformat_target = []
        for target in targets:
            target = re.sub('\_'+self.cytosol, '',target)
            reformat_target.append(target)
        return(reformat_target)

    def get_original_target(self, tmol):
        """Obtain original target"""
        index = None
        for count, target in enumerate(self.targets):
            if tmol in target:
                index =  count
         
        if not index: 
            print ('Could not get target '+tmol)
            return (None)
        else:
            return index
    def calculate_tanimoto_score(self):
        """Calculate similarity for compounds that are not in the database"""
        self.outlier_cpds = set(self.individualtargets) - set(self.all_compounds)
        if self.outlier_cpds:
            self.finaltargets = self.targets
            self.outlier_cpds = self.process_targets(self.outlier_cpds)
            self.retrieve_fingerprints()
            for count ,tmol in enumerate(self.outlier_cpds_fp.keys()):
                if tmol in self.cpd_dict.keys():
                    index = self.get_original_target(tmol+'_'+self.cytosol)
                    if index:
                        new_target = self.extract_db_cpd_ID(tmol)
                        print ("Switching target ID from {} to {}".format(tmol+'_'+self.cytosol, new_target))
                        self.finaltargets[index][0] = new_target
                else:
                    temp = {}
                    for db_cpd in self.db_cpds_fp:
                        score = self.IN.similarity(self.outlier_cpds_fp[tmol], self.db_cpds_fp[db_cpd], 'tanimoto')
                        temp[db_cpd] = score
                    max_score = max(temp.values())
                    max_score_cpd = [temp.keys()[temp.values().index(i)] for i in temp.values() if i == max_score]
                    if max_score == self.threshold_score:
                        new_target = self.extract_db_cpd_ID(max_score_cpd[0])
                        print ("Changed {} to {}".format(tmol+'_'+self.cytosol, new_target))
                        if new_target not in self.individualtargets:
                            index = self.get_original_target(tmol+'_'+self.cytosol)
                            if index:
                                self.finaltargets[index][0] = new_target
                                self.individualtargets[index] = new_target
                        elif new_target in self.individualtargets:
                            print ("New target {} already in list...removing old target {}".format(new_target, tmol+'_'+self.cytosol))
                            index = self.get_original_target(tmol+'_'+self.cytosol)
                            if index:
                                self.finaltargets.remove(self.finaltargets[index])
        else:
            print ("All compounds are in the database")
            self.finaltargets = self.targets

    def retrieve_fingerprints(self):
        """Retrieve fingerprints for database compounds and targets"""
        db_cpds = []
        self.cpd_dict = {}
        for cpd in self.all_compounds:
            originalcpd = deepcopy(cpd)
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
            db_cpds.append(cpd)
            self.cpd_dict.setdefault(cpd, []).append(originalcpd)

        self.outlier_cpds_fp = {}
        for count, tmol in enumerate(self.outlier_cpds):
            try:
                self.outlier_cpds_fp[tmol] = self.INCHI.loadMolecule(tmol).fingerprint('full')
            except indigo.IndigoException:
                print ('Could not get fingerprint for {}'.format(tmol))        
        self.db_cpds_fp = {}
        db_cpds_set = set(db_cpds)
        for db_cpd in db_cpds_set:
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
