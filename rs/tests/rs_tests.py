from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests main rs module'
import re
import os
import unittest
import subprocess

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)
_database = PATH+'/test.db'
_constraints = PATH+'/test.constraints'
_dir = os.path.dirname(__file__)
_executable = os.path.abspath(os.path.join(_dir, os.pardir, 'rs.py'))
_TEST_TARGETS = PATH+'/data2/test_targets1.txt'
_MEDIA = PATH+'/data2/media.txt'

class BRSTests(unittest.TestCase):
    def setUp(self):
        '''Create an instance of main'''
        print ("Initializing tests")

    def tearDown(self):
        '''Delete data structure'''
        print ("Clearing out test suite")

    def test_parse_argument(self):
        '''Tests rs.py'''
        print ("Testing input parameters")
        print ("Testing empty run")

        args = ['python', _executable]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        print ("Testing incompatible modules")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-d_dir', PATH+'/data']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb', _database,
                '-gdbc', _constraints, '-d_dir', PATH+'/data', '-media', _MEDIA]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data', '-ko']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data', '-media', _MEDIA, '-ko']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        print ("Testing compatible modules")
        print ("...Testing generate database")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data', '-rdf', PATH+'/data3']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print ("...Testing generate database with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data4',
                '--inchidb', '-rdf', PATH+'/data3']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()

        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--inchidb', '-rdf', PATH+'/data3']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-d_dir', PATH+'/data4', '--inchidb']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)


        print ("...Testing generate database constraints from preexisting database")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db',
                _database, '-gdbc', _constraints]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print ("...Testing flux balance analysis")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-fba', '-d_dir', PATH+'/data']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-db',
                _database, '-dbc', _constraints, '-fba']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing flux balance analysis with media option")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '-fba', '-media',
                _MEDIA]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing flux balance analysis with reaction knockout option")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db',
                _database, '-dbc', _constraints, '-fba', '-ko']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing flux balance analysis with media, knockout and RDF option")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database, '-dbc', _constraints, '-fba',
                '-media', _MEDIA, '-ko',
                '-rdf', PATH+'/data3']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)
        if os.path.isfile(PATH+'/'+'essentialrxns_output.txt') is True:
            os.remove(PATH+'/essentialrxns_output.txt')
        if os.path.isfile(PATH+'/'+'flux_output.txt') is True:
            os.remove(PATH+'/flux_output.txt')
        if os.path.isfile(PATH+'/'+'fluxKO_output.txt') is True:
            os.remove(PATH+'/fluxKO_output.txt')
        if os.path.isfile(PATH+'/'+'optimal_pathways.txt') is True:
            os.remove(PATH+'/optimal_pathways.txt')
        if os.path.isfile(PATH+'/'+'optimal_pathways.txt') is True:
            os.remove(PATH+'/optimal_pathways.txt')
        if os.path.isfile(PATH+'/'+'theoretical_yield.txt') is True:
            os.remove(PATH+'/theoretical_yield.txt')
        if os.path.isfile(PATH+'/'+'active_metabolism.txt') is True:
            os.remove(PATH+'/active_metabolism.txt')
        if os.path.isfile(PATH+'/'+'flux_individualfluxes_output.txt') is True:
            os.remove(PATH+'/flux_individualfluxes_output.txt')

        if os.path.isfile(PPATH+'/'+'essentialrxns_output.txt') is True:
            os.remove(PPATH+'/essentialrxns_output.txt')
        if os.path.isfile(PPATH+'/'+'fluxKO_output.txt') is True:
            os.remove(PPATH+'/fluxKO_output.txt')
        if os.path.isfile(PPATH+'/'+'flux_output.txt') is True:
            os.remove(PPATH+'/flux_output.txt')
        if os.path.isfile(PPATH+'/'+'optimal_pathways.txt') is True:
            os.remove(PPATH+'/optimal_pathways.txt')
        if os.path.isfile(PPATH+'/'+'optimal_pathways.txt') is True:
            os.remove(PPATH+'/optimal_pathways.txt')
        if os.path.isfile(PPATH+'/'+'theoretical_yield.txt') is True:
            os.remove(PPATH+'/theoretical_yield.txt')
        if os.path.isfile(PPATH+'/'+'active_metabolism.txt') is True:
            os.remove(PPATH+'/active_metabolism.txt')
        if os.path.isfile(PPATH+'/'+'flux_individualfluxes_output.txt') is True:
            os.remove(PPATH+'/flux_individualfluxes_output.txt')
if __name__ == '__main__':
    unittest.main()
