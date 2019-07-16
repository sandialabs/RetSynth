from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Gets IDs from compound chemical formulas or names'

import re
from Pubchem import pubchem_compounds as pc
from copy import deepcopy
class Readfile(object):
    """Reads the input file containing target compounds"""
    def __init__(self, file_name, db, inchidb=False, compartment='cytosol'):
        '''Initialize class'''
        self.file_name = file_name
        self.DB = db
        self.targets = []
        self.inchidb = inchidb
        self.compartment = compartment
        self.PC = pc.PubchemConnector(self.DB)
        self.file_opener()
        self.modelinfo_opener()
        self.get_db_info()

    def file_opener(self):
        '''
        Opens target file and extracts information
        '''
        self.input = {}
        header = {}
        with open(self.file_name) as f:
            line = f.readline()
            if line.startswith('#'):
                line = line.replace("#", "")
                line = line.lower()
                head = line.strip('\n').split('\t')
                for count, h in enumerate(head):
                    header[count] = h.lower()
            else:
                raise ValueError('Input file needs header or # sign in front of header line')

            for line_count, line in enumerate(f):
                self.input[line_count] = {}
                larray = line.strip('\n').split('\t')
                for count, item in enumerate(larray):
                    self.input[line_count][header[count]] = item

    def get_db_info(self):
        '''
        Retrieve database information for the values in the inputfile
        '''
        self.targets = []
        self.ignorerxns = []
        self.includerxns = []
        compartmentID = self.DB.get_compartment(self.compartment)
        if not compartmentID or len(compartmentID) == 0:
            print ('WARNING: No compartment info going to default c0')
            compartmentID = 'c0'
        else:
            compartmentID = compartmentID[0]        
        for count, values in self.input.items():
            temp = []
            if 'compoundid' in values:
                match = re.search('\_\w{1}0$', values['compoundid'])
                if match:
                    temp.extend((values['compoundid'], ''))
                else:
                    temp.extend((values['compoundid']+'_'+compartmentID, ''))                    
            elif 'pubchem' in values:
                cid = values['pubchem']
                db_cpdID = self.PC.get_ID_from_pubchemID(cid, self.inchidb)
                if db_cpdID != '':
                    if isinstance(db_cpdID, list):
                        if len(db_cpdID) == 1:
                            temp.extend((db_cpdID[0], cid))
                        elif len(db_cpdID) > 1:
                            print ('WARNING: Multiple cpd IDs found for one pubchemID automatically choses first ID {}'.format(db_cpdID))
                            temp.extend((db_cpdID[0], cid))
                    else:
                        temp.extend((db_cpdID, cid))
                else:
                    print ('WARNING: no compound ID found for pubchem dentifier {} value not added to list'.format(cid))
            elif 'name' in values:
                name = values['name']
                db_cpdID = self.PC.get_ID_from_name(name)
                if db_cpdID != '':
                    temp.extend((db_cpdID, name))
                else:
                    print ('WARNING: no compound ID found for compound name {} value not added to list'.format(name))
            elif 'inchi' in values:
                temp.extend((values['inchi']+'_'+compartmentID, ''))
            else:
                raise ValueError('Need compound ID, pubchem ID or name of compound in target file')

            if temp:
                if 'organismid' in values.keys():
                    multiplevalues = values['organismid'].split(',')
                    for m in multiplevalues:
                        temp2 = deepcopy(temp)
                        temp2.extend((m, ''))
                        self.targets.append(temp2)
                elif 'organism' in values.keys():
                    multipleorgs = values['organism'].split(',')
                    for org in multipleorgs:
                        if org == 'NA' or org == 'N/A':
                            temp2.extend((org, ''))
                            self.targets.append(temp2)
                        else:
                            try:
                                new_org_format = re.sub(' ', '', org)
                                new_org_format = re.sub('_', '', new_org_format)
                                orgID = self.modelfilename_dict[new_org_format.upper()]
                                temp2 = deepcopy(temp)
                                temp2.extend((orgID, org))
                                self.targets.append(temp2)
                            except KeyError:
                                print ('WARNING: No exact species ID could be found for {} will look to see if multiple strains exist'.format(org))
                                strains = self.strain_check(new_org_format)
                                if strains:
                                    print ('WARNING: Found {} strains of {}, if more than 5 recommend specifying strain in organisms name'.format(len(strains), org))
                                    for orgID in strains:
                                        temp2 = deepcopy(temp)
                                        temp2.extend((orgID, org))
                                        self.targets.append(temp2)
                                else:
                                    print ('WARNING: No species ID could be found for {} will proceed without ID'.format(org))
                                    temp2 = deepcopy(temp)
                                    temp2.extend(('', org))
                                    self.targets.append(temp2)
                	
                else:
                    print ('WARNING: No organism was given or found')
                    temp.extend(('', ''))
                    self.targets.append(temp)

                if 'ignore reactions' in values.keys():
                    multiplereactions = values['ignore reactions'].split(',')
                    self.ignorerxns = self.ignorerxns+multiplereactions

                if 'include reactions' in values.keys():
                    multiplereactions = values['ignore reactions'].split(',')
                    self.includerxns = self.includerxns+multiplereactions

    def strain_check(self,new_org_format):
        '''Extract all strains for a specified organism'''
        strains = []
        for key in self.modelfilename_dict:
            if key.startswith(new_org_format.upper()):
                strains.append(self.modelfilename_dict[key])
        return (strains)

    def modelinfo_opener(self):
        '''
        Get all model information out of the database
        '''
        self.organisms = {}
        modelIDs = self.DB.get_all_fba_models()
        self.modelfilename_dict = {}
        self.modelID_dict = {}
        for file_name in modelIDs:
            fn = re.sub('.xml', '', file_name[1])
            fn = re.sub('_', ' ', fn)
            fn = re.sub('\s+', '', fn)
            self.modelfilename_dict[fn.upper()] = file_name[0]
            self.modelID_dict[file_name[0]] = fn
