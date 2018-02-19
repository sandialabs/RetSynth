from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests main rs module'
import re
import os
import shutil
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
_OUTPUT_DIRECTORY = PATH+'/output_test'
try:
    os.mkdir(_OUTPUT_DIRECTORY)
except OSError:
    shutil.rmtree(_OUTPUT_DIRECTORY)
    os.mkdir(_OUTPUT_DIRECTORY)

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
                _database, '-gdbc', _constraints, '--kbase']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--atlas']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--mine']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--SPRESI']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-k_dir', PATH+'/data']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb', _database,
                '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '-media', _MEDIA]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '-ko']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data',
                 '-media', _MEDIA, '-ko']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)

        print ("Testing compatible modules")
        print ("...Testing generate database")

        ###############TEST GENERATION OF A NOVEL DATABASE###############
        ##TEST GENERATION OF DATABASE NOT USING INCHI IDS##
        print("...Testing kbase database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing kbase & metacyc database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase',
                '--metacyc', '-mc', PATH+'/data7/metabolic-reactions.xml',
                 '-tf', PATH+'/data7/MetaCyc.aliases',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing kegg (only) database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--kegg', '-keggnunorganisms', '1',
                '-keggnunorganismpaths', '1',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing kegg with kbase database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--kegg', '-keggnunorganisms', '1',
                '-keggnunorganismpaths', '1', '--kbase', '-k_dir', PATH+'/data',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing spresi database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--SPRESI',
                '-s_dir', PATH+'/data3', '--kbase', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing atlas database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '--atlas',
                '-a_dir', PATH+'/data5', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing mine database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase',
                '--mine', '-m_dir', PATH+'/data6', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing spresi & mine database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data',
                 '--SPRESI', '-s_dir', PATH+'/data3', '--kbase', '--mine',
                 '-m_dir', PATH+'/data6', '-op', _OUTPUT_DIRECTORY]

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing spresi, mine & atlas database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--SPRESI',
                '-s_dir', PATH+'/data3', '--kbase', '--mine', '-m_dir', PATH+'/data6',
                '--atlas', '-a_dir', PATH+'/data5', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing metacyc, spresi, mine & atlas database generation without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--SPRESI',
                '-s_dir', PATH+'/data3', '--kbase', '--mine', '-m_dir', PATH+'/data6',
                '--atlas', '-a_dir', PATH+'/data5', '--metacyc', '-mc',
                  PATH+'/data7/metabolic-reactions.xml', '-tf', PATH+'/data7/MetaCyc.aliases',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        ##TEST GENERATION OF DATABASE USING INCHI IDS##
        print("...Testing kbase database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase',
                '--inchidb','-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing kbase & metacyc database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase',
                '--metacyc', '-mc', PATH+'/data7/metabolic-reactions.xml',
                  '-tf', PATH+'/data7/MetaCyc.aliases', '--inchidb', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)
        print("...Testing kegg (only) database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--kegg', '-keggnunorganisms', '1',
                '-keggnunorganismpaths', '1', '--inchidb',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing kegg with kbase database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '--kegg', '-keggnunorganisms', '1',
                '-keggnunorganismpaths', '1', '--kbase', '-k_dir', PATH+'/data',
                '--inchidb', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing spresi database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '--SPRESI',
                '-s_dir', PATH+'/data3', '-m_dir', PATH+'/data6', '--inchidb',
                 '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing atlas database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '--atlas',
                '-a_dir', PATH+'/data5', '--inchidb','-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing mine database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '--mine',
                '-m_dir', PATH+'/data6', '--inchidb', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing spresi & mine database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '--SPRESI',
                '-s_dir', PATH+'/data3', '-m_dir', PATH+'/data6', '--mine',
                '-m_dir', PATH+'/data6', '--inchidb', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print("...Testing spresi, mine & atlas database generation with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data', '--kbase', '--SPRESI',
                '-s_dir', PATH+'/data3', '-m_dir', PATH+'/data6', '--mine',
                '-m_dir', PATH+'/data6', '--atlas',
                '-a_dir', PATH+'/data5', '--inchidb', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        ##############TEST ADDITION OF A NEW DATABASE TO A PRE-EXISTING DATABASE###############

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data4',
                '--kbase', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()

        print ("...Testing metacyc addition to database without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--metacyc', '-mc', PATH+'/data7/metabolic-reactions.xml',
                '-tf', PATH+'/data7/MetaCyc.aliases', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing kegg addition to database without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--kegg', '-keggnunorganisms', '1',
                '-keggnunorganismpaths', '1', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing spresi addition to database without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--SPRESI', '-s_dir', PATH+'/data3', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing mine addition to database without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--mine', '-m_dir', PATH+'/data6', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing atlas addition to database without inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--atlas', '-a_dir', PATH+'/data5', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-k_dir', PATH+'/data4',
                '--inchidb', '--kbase', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()

        print ("...Testing metacyc addition to database with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--inchidb',
                '--metacyc', '-mc', PATH+'/data7/metabolic-reactions.xml',
                '-tf', PATH+'/data7/MetaCyc.aliases',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing kegg addition to database with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--kegg', '-keggnunorganisms', '1',
                '-keggnunorganismpaths', '1', '--inchidb', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)


        print ("...Testing spresi ddition to database with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--inchidb', '--SPRESI', '-s_dir', PATH+'/data3',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing mine ddition to database with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--inchidb', '--mine', '-m_dir', PATH+'/data6',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing atlas ddition to database with inchi")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '--inchidb', '--atlas', '-a_dir', PATH+'/data5',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-op', _OUTPUT_DIRECTORY, '--kbase',
                '-k_dir', PATH+'/data4', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print ("...Testing generate database constraints from preexisting database")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db',
                _database, '-gdbc', _constraints, '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)

        print ("...Testing flux balance analysis")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-gdb',
                _database, '-gdbc', _constraints, '-fba', '-k_dir', PATH+'/data', '--kbase',
                '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        args = ['python', _executable, '-t', _TEST_TARGETS, '-db',
                _database, '-dbc', _constraints, '-fba', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing flux balance analysis with media option")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '-fba', '-media',
                _MEDIA, '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing flux balance analysis with reaction knockout option")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db',
                _database, '-dbc', _constraints, '-fba', '-ko', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

        print ("...Testing flux balance analysis with media, knockout and RDF option")
        args = ['python', _executable, '-t', _TEST_TARGETS, '-db', _database,
                '-dbc', _constraints, '-fba', '-media', _MEDIA, '-ko', '--SPRESI',
                '-s_dir', PATH+'/data3', '-op', _OUTPUT_DIRECTORY]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        if os.path.isfile(_database) is True:
            os.remove(_database)
        if os.path.isfile(_constraints) is True:
            os.remove(_constraints)
        try:
            shutil.rmtree(_OUTPUT_DIRECTORY)
        except:
            pass

if __name__ == '__main__':
    unittest.main()
