from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'builds tables MINE database from msp files'
import os
import re
import glob
import sqlite3
from sys import platform
from tqdm import tqdm
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_inchi
elif platform == "win32" or platform == "win64":
    from indigopython130_win import indigo
    from indigopython130_win import indigo_inchi
from Database import query as Q

class BuildMINEdb(object):
    '''Adds or builds MINE database to metabolic database'''
    def __init__(self, dumpdirectory, database, inchidb, rxntype):
        '''initialize'''
        self.database = database
        self.rxntype = rxntype
        self.inchidb = inchidb
        self.IN = indigo.Indigo()
        self.INCHI = indigo_inchi.IndigoInchi(self.IN)
        self.compound_dict = {}
        self.reaction_dict = {}
        self.compound_dict_temp = {}
        self.storerxns = set()
        self.count_reaction = 0
        msp_files = glob.glob(os.path.join(dumpdirectory, '*'))
        for filename in msp_files:
            self.open_mspfile(filename)
            self.generate_reactions()
        self.fill_database()

    def extract_cpd_information(self, compoundid, INFO, tp):
        '''Get compound information'''
        if tp.startswith(INFO):
            tparray = tp.split()
            if INFO not in self.compound_dict[compoundid]:
                self.compound_dict[compoundid][INFO] = tparray[1]
                if INFO == 'SMILES' and self.inchidb is True:
                    try:
                        mol = self.IN.loadMolecule(tparray[1])
                        if mol:
                            inchi = self.INCHI.getInchi(mol)
                            # fp = mol.fingerprint('full')
                            # buffer = fp.toBuffer()
                            # buffer_array = [str(i) for i in buffer]
                            # buffer_string = ','.join(buffer_array)
                            cf = mol.grossFormula()
                            cf = re.sub(' ', '', cf)                           
                            if inchi:
                                self.compound_dict[compoundid]['INCHI'] = inchi
                                # self.compound_dict[compoundid]['FP'] = buffer_string
                                self.compound_dict[compoundid]['CF'] = cf
                        else:
                            pass
                    except indigo.IndigoException:
                        pass
    def extract_source_information(self, compoundid, tp):
        '''extract operator information (EC number)'''
        if tp.startswith('Sources'):
            tp = re.sub('Sources: ', '', tp)
            tp = re.sub('\[|\]', '', tp)
            operators = tp.split('}, {')
            for operator in operators:
                operator = re.sub('\{|\}', '', operator)
                operator = re.sub('u\'', '', operator)
                operator = re.sub('\'', '', operator)
                operator = re.sub('Operators: ', '', operator)
                #operator = re.sub('Compound: ', '', operator)
                ops1 = operator.split(', Compound: ')
                enzymes = ops1[0].split(', ')
                try:
                    compound = ops1[1]
                    for enzyme in enzymes:
                        self.compound_dict[compoundid]['ENZYME'].setdefault(enzyme,
                                                                            set()).add(compound)
                except IndexError:
                    print ('Index error for '+str(compoundid)+' '+str(operator)+' '+str(ops1))

    def add2dictionary(self, temp):
        '''Adds compound to file dictionary'''
        compoundid = ''
        for tp in temp:
            if tp.startswith('_id'):
                tparray = tp.split()
                if tparray[1] not in self.compound_dict:
                    self.compound_dict[tparray[1]] = {}
                compoundid = tparray[1]
        if 'ENZYME' not in self.compound_dict[compoundid]:
            self.compound_dict[compoundid]['ENZYME'] = {}

        for tp in temp:
            self.extract_cpd_information(compoundid, 'SMILES', tp)
            self.extract_cpd_information(compoundid, 'Name', tp)
            self.extract_cpd_information(compoundid, 'Inchikey', tp)
            self.extract_cpd_information(compoundid, 'MINE_id', tp)
            self.extract_source_information(compoundid, tp)
        return []

    def open_mspfile(self, filename):
        '''Opens and reads msp files'''
        print ('STATUS: Reading '+filename)
        numlines = len(open(filename).readlines())
        with open(filename) as fin:
            cpd = []
            for count_line, line in enumerate(tqdm(fin)):
                count_line = count_line+1
                line = line.strip()
                if line != '':
                    cpd.append(line)
                else:
                    cpd = self.add2dictionary(cpd)
                if int(count_line) is int(numlines):
                    cpd = self.add2dictionary(cpd)

    def fill_reaction_components_dict(self, rxn, compound, typecpd):
        '''Fill reaction info (substrates) in to dictionary'''
        if self.inchidb:
            try:
                value = self.compound_dict[compound]
                try:
                    self.reaction_dict[rxn][typecpd].append(self.compound_dict[compound]['INCHI'])
                except KeyError:
                    self.reaction_dict[rxn][typecpd].append(compound)
            except KeyError:
                self.compound_dict_temp[compound] = {}
                self.compound_dict_temp[compound]['Name'] = 'None'
                self.compound_dict_temp[compound]['SMILES'] = 'None'
                self.compound_dict_temp[compound]['Inchikey'] = 'None'
                self.compound_dict_temp[compound]['MINE_id'] = 'None'
                self.compound_dict_temp[compound]['ENZYME'] = {}
                self.reaction_dict[rxn][typecpd].append(compound)
        else:
            try:
                value = self.compound_dict[compound]
                self.reaction_dict[rxn][typecpd].append(compound)
            except KeyError:
                self.compound_dict_temp[compound] = {}
                self.compound_dict_temp[compound]['Name'] = 'None'
                self.compound_dict_temp[compound]['SMILES'] = 'None'
                self.compound_dict_temp[compound]['Inchikey'] = 'None'
                self.compound_dict_temp[compound]['MINE_id'] = 'None'
                self.compound_dict_temp[compound]['ENZYME'] = {}
                self.reaction_dict[rxn][typecpd].append(compound)
    
    def generate_reactions(self):
        '''Get reactions from MINE files'''
        print ('STATUS identify reactions in file')
        for compound in self.compound_dict:
            for EC in self.compound_dict[compound]['ENZYME']:
                reactants = self.compound_dict[compound]['ENZYME'][EC]
                full_rxn = ', '.join(reactants)+'\t'+EC+'\t'+compound
                if full_rxn not in self.storerxns:
                    self.storerxns.add(full_rxn)
                    self.count_reaction += 1
                    rxn = 'rxn'+str(self.count_reaction)+'_m'
                    self.reaction_dict[rxn] = {}
                    self.reaction_dict[rxn]['ENZYME'] = EC
                    self.reaction_dict[rxn]['reactants'] = []
                    self.reaction_dict[rxn]['products'] = []
                    self.fill_reaction_components_dict(rxn, compound, 'products')
                    for reactant in reactants:
                        self.fill_reaction_components_dict(rxn, reactant, 'reactants')
                else:
                    print (full_rxn+'  already present therefore skipping')
        self.compound_dict.update(self.compound_dict_temp)

    def fill_database(self):
        '''Generate arrays of database information and
        fill database with information'''
        reactions = []
        reaction_compound = []
        reaction_gene = []
        reaction_protein = []
        model_reaction = []
        reaction_reversibility = []
        original_cpds = []
        compound = []
        model_compound = []
        print ('STATUS: Load into database')
        cnx = sqlite3.connect(self.database)
        # conn.text_factory = str
        # cnx = conn.cursor()        
        cnx.execute("PRAGMA synchronous = OFF")
        cnx.execute("PRAGMA journal_mode = OFF")
        cnx.commit()
        cytosol = '_c0'
        QC = cnx.execute("""SELECT ID from model where ID=?""", ('MINE',))
        results = QC.fetchall()
        if not results:
            cnx.execute("INSERT INTO model VALUES (?,?)", ('MINE', 'MINE_database'))
            QC = cnx.execute('''SELECT DISTINCT cluster_num FROM cluster''')
            hits = QC.fetchall()
            uniq_clusters = [i[0] for i in hits]
            cnx.execute("INSERT INTO cluster VALUES (?,?)", (len(uniq_clusters)+1, 'MINE'))
            cnx.commit()

        for reaction in self.reaction_dict:
            for reactant in self.reaction_dict[reaction]['reactants']:
                reaction_compound.append((reaction+cytosol, reactant+cytosol, 0, 1, 0))
            for product in self.reaction_dict[reaction]['products']:
                reaction_compound.append((reaction+cytosol, product+cytosol, 1, 1, 0))
            reaction_reversibility.append((reaction+cytosol, 'false'))
            reactions.append((reaction+cytosol, 'None', 'None', self.rxntype))
            reaction_protein.append((reaction+cytosol, 'MINE',
                                     self.reaction_dict[reaction]['ENZYME']))
            reaction_gene.append((reaction+cytosol, 'MINE', 'None'))
            model_reaction.append((reaction+cytosol, 'MINE', 'false'))
        cnx.executemany("INSERT INTO reaction_reversibility VALUES (?,?)", reaction_reversibility)
        cnx.executemany("INSERT INTO reaction VALUES (?,?,?,?)", reactions)
        cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", reaction_protein)
        cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", reaction_gene)
        cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", model_reaction)
        cnx.executemany("INSERT INTO reaction_compound VALUES (?,?,?,?,?)", reaction_compound)

        cnx.commit()
        DB = Q.Connector(self.database)
        original_cpds_in_db = DB.get_all_compounds()

        temp_inchi = set()
        for cpd in self.compound_dict:
            try:
                cpd_inchi = self.compound_dict[cpd]['INCHI']+cytosol
                if cpd_inchi not in original_cpds_in_db:
                    if self.inchidb:
                        original_cpds.append((cpd+cytosol,
                                              self.compound_dict[cpd]['INCHI']+cytosol))
                    if self.compound_dict[cpd]['INCHI']+cytosol  not in temp_inchi:
                        model_compound.append((self.compound_dict[cpd]['INCHI']+cytosol, 'MINE'))
                        try:
                            compound.append((self.compound_dict[cpd]['INCHI']+cytosol,
                                             self.compound_dict[cpd]['Name'], 'c0', 'None',
                                             self.compound_dict[cpd].get('CF', 'None'),
                                             'None'))
                        except KeyError:
                            compound.append((self.compound_dict[cpd]['INCHI']+cytosol,
                                             'None', 'c0', 'None',
                                              self.compound_dict[cpd].get('CF', 'None'), 
                                              'None'))
                        temp_inchi.add(cpd_inchi)
            except KeyError:
                model_compound.append((cpd+cytosol, 'MINE'))
                try:
                    compound.append((cpd+cytosol, self.compound_dict[cpd]['Name'], 'c0', 'None', 'None', 'None'))
                except KeyError:
                    compound.append((cpd+cytosol, 'None', 'c0', 'None', 'None', 'None'))
        if self.inchidb:
            cnx.executemany("INSERT INTO original_db_cpdIDs VALUES (?,?)", original_cpds)
        cnx.executemany("INSERT INTO model_compound VALUES (?,?)", model_compound)
        cnx.executemany("""INSERT INTO compound VALUES (?,?,?,?,?,?)""", compound)
        cnx.commit()
