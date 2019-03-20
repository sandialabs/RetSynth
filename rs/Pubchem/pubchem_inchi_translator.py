from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Gets InChis for compounds in database'

import re
import httplib
import urllib2
import pubchempy as pcp

class CompoundTranslator(object):
    """ Converts compound IDs to their InChi"""

    def translate(self, compound_name):
        '''
        Retrieve InChi\'s for compounds
        '''
        cas_value = 'None'
        compound_name = re.sub('\_\w{1}0$', '', compound_name)
        compound_name = re.sub('_', ' ', compound_name)
        self.get_inchi(compound_name)

        if len(self.IDs) == 0:
            compound_name = re.sub(' ', '-', compound_name)
            self.get_inchi(compound_name)

        if len(self.IDs) == 0:
            compound_name = compound_name+'+'
            self.get_inchi(compound_name)

        if len(self.IDs) == 0:
            compound_name = compound_name+'-'
            self.get_inchi(compound_name)

        if len(self.IDs) > 0:
            allsynomyms = self.IDs[0].synonyms
            for syn in allsynomyms:
                if syn.startswith('CAS'):
                    cas_value = re.sub('CAS-', '', syn)
            return(self.IDs[0].inchi, self.IDs[0].iupac_name, cas_value)
        else:
            return(None, None, cas_value)

    def get_inchi(self, compound_name):
        '''Attempt to get inchi for a compound'''
        try:
            self.IDs = pcp.get_compounds(compound_name, 'name')
        except (pcp.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError, ValueError):
            self.IDs = []
            print ('WARNING: could not get info for {}...Errored out'.format(compound_name))        
