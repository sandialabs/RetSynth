from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'read file containing starting compounds'

def readfile_startcompounds(file_name, compartment='c0'):
    """Reads the input file containing starting compounds"""
    start_cpds = []
    with open(file_name) as f:
        line = f.readline()
        if line.startswith('#'):
            line = line.replace("#", "")
            line = line.lower()
            head = line.strip('\n').split('\t')

        for line in f:
            line = line.strip('\n')
            if line.endswith(compartment):
                pass
            else:
                line = line+'_'+compartment
            start_cpds.append(line)
    return start_cpds
