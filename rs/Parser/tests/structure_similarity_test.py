from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests structure_similarity code'

import os
import unittest
from Parser import structure_similarity as ss

class StructureSimilarityTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_structure_similarity(self):
        database_cpds = []
        targets = []
        with open ('database_compounds.txt') as fin:
            for line in fin:
                line = line.strip()
                database_cpds.append(line)

        with open('testcompounds.txt') as fin:
            for line in fin:
                line = line.strip()
                if line.startswith('#'):
                    pass
                else:
                    targets.append([line, '','',''])

        SIM = ss.TanimotoStructureSimilarity(targets, database_cpds,'c0', 'e0')
        self.assertIn(['InChI=1S/C5H12O/c1-3-4-5(2)6/h5-6H,3-4H2,1-2H3/t5-/m1/s1_c0', '','',''], SIM.finaltargets)
        self.assertNotIn(['InChI=1S/C5H12O/c1-3-4-5(2)6/h5-6H,3-4H2,1-2H3_c0', '','',''], SIM.finaltargets)
        self.assertIn(['InChI=1S/C14H28O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14(15)16/h2-13H2,1H3,(H,15,16)/p-1_p0', '','',''], SIM.finaltargets)
        self.assertNotIn(['InChI=1S/C14H28O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14(15)16/h2-13H2,1H3,(H,15,16)/p-1_c0', '','',''], SIM.finaltargets)
        self.assertIn(['InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3_c0', '','',''], SIM.finaltargets)
        self.assertIn(['InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3_c0', '','',''], SIM.finaltargets)
        self.assertEqual(len(SIM.finaltargets), 39)

if __name__ == '__main__':
    unittest.main()

