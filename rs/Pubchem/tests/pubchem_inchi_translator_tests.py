from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests pubchem compounds inchi collection'

import unittest
from Pubchem import pubchem_inchi_translator as pit

class Tests_pubchem_inchi_translator(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")
        self.CT = pit.CompoundTranslator()

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
        del self.CT

    def test_pit(self):
        """Test inchi translator"""
        print ("...Testing inchi translator")
        inchi2, iupac_name = self.CT.translate('H_c0')
        self.assertEqual('InChI=1S/p+1', inchi2)
        self.assertEqual('hydron', iupac_name)

        inchi2, iupac_name = self.CT.translate('H2O_c0')
        self.assertEqual('InChI=1S/H2O/h1H2', inchi2)
        self.assertEqual('oxidane', iupac_name)

        inchi2, iupac_name = self.CT.translate('L_Valine_e0')
        self.assertEqual('InChI=1S/C5H11NO2/c1-3(2)4(6)5(7)8/h3-4H,6H2,1-2H3,(H,7,8)/t4-/m0/s1',
                         inchi2)
        self.assertEqual('(2S)-2-amino-3-methylbutanoic acid', iupac_name)

if __name__ == '__main__':
    unittest.main()
