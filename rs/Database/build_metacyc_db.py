from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Translates metacyc into database'

import os
import re
import sqlite3
import urllib2
import httplib
import pubchempy
from tqdm import tqdm
from Database import query as Q
from bs4 import BeautifulSoup, SoupStrainer
import lxml
from Pubchem import pubchem_inchi_translator as pit
from sys import platform
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_inchi
elif platform == "win32" or platform == "win64":
    from indigopython130_win import indigo
    from indigopython130_win import indigo_inchi

PATH = os.path.dirname(os.path.abspath(__file__))
KEGG = 'http://rest.kegg.jp/'

def verbose_print(verbose, line):
    if verbose:
        print(line)

def extract_KEGG_data(url):
    '''Extract Kegg db info'''
    print (url)
    try:
        data = urllib2.urlopen(url).read()
        darray = data.split('\n')
        return darray
    except (httplib.BadStatusLine, urllib2.URLError):
        return None

def get_inchi_from_kegg_ID(cpd):
    darray = extract_KEGG_data(KEGG+'get/'+cpd)
    inchicpd = None
    cas = None
    if darray:
        for value in darray:
            array = value.split()
            if 'PubChem:' in array:
                index = array.index('PubChem:')
                sid = array[index+1]
                try:
                    substance = pubchempy.Substance.from_sid(sid)
                    substance_cids = substance.cids
                    if substance_cids:
                        try:
                            compounds = pubchempy.get_compounds(substance_cids[0])
                            if compounds:
                                inchicpd = compounds[0].inchi
                        except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                            pass
                except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                    print ('WARNING: Could not get substance for {} {}'.format(sid, cpd))
                    pass
            if 'CAS:' in array:
                index = array.index('CAS:')
                cas = array[index+1]
    return (inchicpd, cas)

class Translate(object):
    """Translates metacyc compound and reaction IDs to Kbase compound and reaction IDs
    if a kbase database is used (ONLY WORKS WITH KBASE, ModelSeed/patric)"""
    def __init__(self, DBPath, file_name, inchidb, rxntype, verbose, add=True):
        # print (translation_file)
        # print (file_name)
        self.file_name = file_name
        self.inchidb = inchidb
        self.cnx = sqlite3.connect(DBPath)
        self.cnx.execute("PRAGMA synchronous = OFF")
        self.cnx.execute("PRAGMA journal_mode = OFF")
        self.DB = Q.Connector(DBPath)
        self.rxntype = rxntype
        self.verbose = verbose
        self.BIOCYC_translator = {}
        self.BIOCYC_translator['rxn'] = {}
        self.BIOCYC_translator['compound'] = {}

        with open(PATH+'/data/MetaCyc.aliases') as file_name:
            line = file_name.readline()
            if line.startswith('#'):
                pass
            for line in file_name:
                larray = line.strip('\n').split('\t')
                larray[0] = re.sub('\.\w+$', '', larray[0])
                if larray[1]:
                    self.BIOCYC_translator['rxn'].setdefault(larray[0], []).append(larray[1])
                elif larray[2]:
                    self.BIOCYC_translator['compound'].setdefault(larray[0], []).append(larray[2])

        verbose_print(self, "STATUS: {} BIOCYC IDS".format(len(self.BIOCYC_translator['rxn'])+len(self.BIOCYC_translator['compound'])))
        if add:
            self.add_metacyc_to_db()

    def add_metacyc_to_db(self):
        '''
        Adds metacyc information to the database
        '''
        MC = MetaCyc(self.DB, self.inchidb, self.cnx, self.verbose)
        MC.read_metacyc_file(self.BIOCYC_translator, self.file_name)

        def add_new_cluster_info(MC_mi):
            '''
            Adds updated metacyc cluster info to database
            '''
            Q = self.cnx.execute('''INSERT INTO model VALUES (?,?)''', (MC_mi, self.file_name))
            Q = self.cnx.execute('''SELECT DISTINCT cluster_num FROM cluster''')
            hits = Q.fetchall()
            uniq_clusters = [i[0] for i in hits]
            self.cnx.execute("INSERT INTO cluster VALUES (?,?)", (len(uniq_clusters)+1, MC_mi))
            self.cnx.commit()

        def add_new_cpd_rxn_gene_info(MC):
            '''
            Adds metacyc compound, reaction, gene and protein info to database
            '''
            self.cnx.executemany("INSERT INTO model_compound VALUES (?,?)", MC.modelcompounds)
            self.cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", MC.modelreactions)
            self.cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", MC.genelist)
            self.cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", MC.proteinlist)
            self.cnx.commit()

        def add_cpds_2_db(allcompounds):
            '''add compounds to database'''
            count_compound = 0
            for cpd in allcompounds:
                Q = self.cnx.execute("SELECT ID FROM compound WHERE ID = ?", (cpd[0],))
                result = Q.fetchone()
                if result is None:
                    count_compound += 1
                    self.cnx.execute("""INSERT INTO compound VALUES (?,?,?,?,?,?)""",
                                     (cpd[0], cpd[1], cpd[2], cpd[3], cpd[4], cpd[5]))
            self.cnx.commit()
            print ('STATUS: {} added MetaCyc compounds'.format(count_compound))

        def add_rxns_2_db(allreactions):
            '''add reaciton to database'''
            count_reactions = 0
            count_reactions_compounds = 0

            for rxn in allreactions:
                Q = self.cnx.execute("SELECT ID FROM reaction WHERE ID = ?", (rxn[0],))
                result = Q.fetchone()
                if result is None:
                    count_reactions += 1
                    count_reactions_compounds += 1
                    self.cnx.execute("INSERT INTO reaction VALUES (?,?,?,?)", (rxn[0], rxn[1], rxn[2], self.rxntype))
                    # self.cnx.commit()
                    for rxn_compound in MC.all_reaction_compound[rxn[0]]:
                        self.cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                         (rxn_compound[0], rxn_compound[1],
                                          rxn_compound[2], rxn_compound[3], 0))
            self.cnx.commit()
            print ('STATUS: {} added MetaCyc reactions'.format(count_reactions))
            print ('STATUS: {} added MetaCyc reaction_compounds'.format(count_reactions_compounds))

        def determine_rxn_reversibility(reactionreversibility):
            '''set reversibility'''
            count_reversibility = 0
            count_reversibility_add = 0
            for rxn_revers in reactionreversibility:
                Q = self.cnx.execute("""SELECT reaction_ID FROM reaction_reversibility \
                                    WHERE reaction_ID=?""", (rxn_revers[0],))
                result = Q.fetchone()
                if result is None:
                    count_reversibility_add += 1
                    self.cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)",
                                     (rxn_revers[0], rxn_revers[1]))
                    # self.cnx.commit()
                else:
                    Q = self.cnx.execute("""SELECT is_reversible FROM reaction_reversibility \
                                         WHERE reaction_ID = ?""", (rxn_revers[0],))
                    result = Q.fetchone()
                    if result[0] != rxn_revers[1] and result[0] == 'false':
                        count_reversibility += 1
                        command = '''UPDATE reaction_reversibility SET is_reversible = 'true' \
                                   WHERE reaction_ID = ?'''
                        self.cnx.execute(command, (rxn_revers[0],))
            self.cnx.commit()
            print ('STATUS: {} changed reversibility'.format(count_reversibility))
            print ('STATUS: {} added MetaCyc reversibility'.format(count_reversibility_add))
   
        def compile_all_cpd_rxn_info(MC):
            '''
            Compiles and adds new compound and reaction information to tables in databases
            that have all reactions and compounds, no repeats of information allowed in these
            tables
            '''
            allreactions = []
            allcompounds = []
            reactionreversibility = []

            for key, value in MC.all_rxnreversibility.iteritems():
                reactionreversibility.append((key, value))
            for key, value in MC.all_rxns.iteritems():
                allreactions.append((key, value[0], value[1]))
            for key, value in MC.all_cpds.iteritems():
                allcompounds.append((key, value[0], value[1], value[2], value[3], value[4]))

            add_cpds_2_db(allcompounds)
            add_rxns_2_db(allreactions)
            determine_rxn_reversibility(reactionreversibility)

        add_new_cluster_info(MC.mi)
        add_new_cpd_rxn_gene_info(MC)
        compile_all_cpd_rxn_info(MC)

#######################READ IN METACYC DATABASE######################
class MetaCyc(object):
    """Opens and parses metacyc xml file"""
    def __init__(self, DB, inchidb, cnx, verbose):
        '''initialize'''
        print ('STATUS loading in metacyc file')
        self.DB = DB
        self.inchidb = inchidb
        self.cnx = cnx
        self.verbose = verbose
        if self.inchidb:
            self.IN = indigo.Indigo()
            self.INCHI = indigo_inchi.IndigoInchi(self.IN)
            self.prom_cpds = self.get_promiscuous_cpds(PATH+'/data/promiscuous_cpds_inchi.txt')
        else:
            self.prom_cpds = self.get_promiscuous_cpds(PATH+'/data/promiscuous_cpds_kegg.txt')

    def get_promiscuous_cpds(self, file_name):
        prom_cpds=set()
        with open(file_name) as fin:
            for line in fin:
                line = line.strip()
                prom_cpds.add(line)
        return prom_cpds

    def compound_translator(self, compound_ID, biocyc_ID, inchi_ID, KEGG_ID, name, compartment):
        '''
        Checks if metacyc compound is in database, if it is not it adds it
        '''
        if KEGG_ID:
            self.fill_compound_arrays(compound_ID, inchi_ID, KEGG_ID, KEGG_ID, name, compartment)
        elif biocyc_ID:
            key = re.sub('\|', '', biocyc_ID)
            try:
                cpdID_list = self.BIOCYC_translator['compound'][key]
                for cpdID in cpdID_list:
                    cpdIDs = cpdID.split('|')
                    if len(cpdIDs) > 1:
                        temp_cpdIDs = []
                        for c in cpdIDs:
                            if c+'_'+compartment in self.all_compounds:
                                temp_cpdIDs.append(c)
                        if len(temp_cpdIDs) == 1:
                            cpdID = temp_cpdIDs[0]
                        elif len(temp_cpdIDs) == len(cpdIDs):
                            cpdID = temp_cpdIDs[0]
                        elif len(temp_cpdIDs) == 0:
                            cpdID = cpdIDs[0]
                        else:
                            print ('WARNING: compound ID {} issue'.format(cpdID))
                    self.fill_compound_arrays(compound_ID, inchi_ID, cpdID, str(KEGG_ID), name, compartment)
            except KeyError:
                self.fill_compound_arrays(compound_ID, inchi_ID, compound_ID, str(KEGG_ID), name, compartment)
        else:
            self.fill_compound_arrays(compound_ID, inchi_ID, compound_ID, str(KEGG_ID), name, compartment)

    def get_fp_cf_info(self, inchi):
        inchi = re.sub('_'+'c0'+'$', '', inchi)
        inchi = re.sub('_'+'e0'+'$', '', inchi)
        try:
            mol = self.INCHI.loadMolecule(inchi)
            # fp = mol.fingerprint('full')
            # buffer = fp.toBuffer()
            # buffer_array = [str(i) for i in buffer]
            # buffer_string = ','.join(buffer_array)
            cf = mol.grossFormula()
            cf = re.sub(' ', '', cf)
        except indigo.IndigoException:
            print ('STATUS: Could not load inchi molecule {}'.format(inchi))
            cf = 'None'
        return (cf)

    def fill_compound_arrays(self, ID, inchi, cpdID, KEGG_ID, name, compartment):
        '''
        Determines whether or not to add compound to compound arrays which will
        later be added to the sqlite database
        '''
        if not self.inchidb:
            self.arrays(cpdID+'_'+compartment, compartment, name, KEGG_ID, 'None', 'None')
            self.store_compounds[ID] = cpdID+'_'+compartment
        else:
            if KEGG_ID != 'None' or KEGG_ID is None:
                cpd_inchi, cas = get_inchi_from_kegg_ID(KEGG_ID)
            else:
                cas=None
            if not inchi:
                cpdID = self.check_db(cpdID+'_'+compartment, KEGG_ID)
                if cpdID.startswith('InChI'):
                    cpdID_temp = re.sub('_'+compartment+'$', '', cpdID)
                    cf = self.get_fp_cf_info(cpdID_temp)
                    self.arrays(cpdID, compartment, name, KEGG_ID, cf, str(cas))
                    self.store_compounds[ID] = cpdID
                else:
                    if KEGG_ID != 'None' or KEGG_ID is None:
                        count = 0
                        while count < 3:
                            cpd, cas = get_inchi_from_kegg_ID(KEGG_ID)
                            if cpd:
                                count = 4
                            else:
                                count+=1
                                print ('WARNING: No INCHI was found from {} will try again for {} time'.format(KEGG_ID, count))                        

                        if cpd:
                            cf = self.get_fp_cf_info(cpd)
                            self.arrays(cpd+'_'+compartment, compartment, name, KEGG_ID, cf, str(cas))
                            self.store_compounds[ID] = cpd+'_'+compartment
                        else:    
                            inchi2, iupac_name, cas = self.CT.translate(name)
                            if not inchi2:
                                self.arrays(cpdID, compartment, name, KEGG_ID, 'None', str(cas))
                                self.store_compounds[ID] = cpdID
                            else:
                                cf = self.get_fp_cf_info(inchi2)
                                self.arrays(inchi2+'_'+compartment, compartment, name, KEGG_ID, cf,  str(cas))
                                self.store_compounds[ID] = inchi2+'_'+compartment
                    else:
                        inchi2, iupac_name, cas = self.CT.translate(name)
                        if not inchi2:
                            self.arrays(cpdID, compartment, name, KEGG_ID, 'None',  str(cas))
                            self.store_compounds[ID] = cpdID
                        else:
                            cf = self.get_fp_cf_info(inchi2)
                            self.arrays(inchi2+'_'+compartment, compartment, name, KEGG_ID, cf,  str(cas))
                            self.store_compounds[ID] = inchi2+'_'+compartment

            else:
                cpdID = self.check_db(cpdID+'_'+compartment, KEGG_ID)
                if cpdID.startswith('InChI'):
                    cpdID_temp = re.sub('_'+compartment+'$', '', cpdID)
                    cf = self.get_fp_cf_info(cpdID_temp)
                    self.arrays(cpdID, compartment, name, KEGG_ID, cf, str(cas))
                    self.store_compounds[ID] = cpdID
                else:
                    try:
                        QC = self.cnx.execute('''SELECT * from compound where ID=?''', (cpdID,))
                        try:
                            result = QC.fetchone()[0]
                            cf = self.get_fp_cf_info(inchi)
                            self.arrays(cpdID, compartment, name, KEGG_ID, cf, result[-1])
                            self.store_compounds[ID] = cpdID
                        except TypeError:
                            self.cnx.execute("INSERT INTO original_db_cpdIDs VALUES (?,?)",
                                             (cpdID, inchi+'_'+compartment))
                            cf = self.get_fp_cf_info(inchi)
                            self.arrays(inchi+'_'+compartment, compartment, name, KEGG_ID, cf, str(cas))
                            self.store_compounds[ID] = inchi+'_'+compartment
                    except sqlite3.OperationalError:
                        # self.cnx.execute("INSERT INTO original_db_cpdIDs VALUES (?,?)",
                        #                  (cpdID+'_'+compartment, inchi+'_'+compartment))
                        cf = self.get_fp_cf_info(inchi)
                        self.arrays(inchi+'_'+compartment, compartment, name, KEGG_ID, cf, str(cas))
                        self.store_compounds[ID] = inchi+'_'+compartment

    def arrays(self, cpdID, compartment, name, KEGG_ID, chemicalformula, cas):
        '''
        Adds compound to compound arrays
        '''
        if cpdID not in self.all_cpds:
            self.all_cpds[cpdID] = [name, compartment, KEGG_ID, chemicalformula, cas]
        if (cpdID, self.mi) not in self.modelcompounds:
            self.modelcompounds.append((cpdID, self.mi))

    def rxn_translator(self, reaction_ID, temp_all_rxn_compound, revers, name, genes, proteins, kegg):
        '''
        Checks if metacyc reaction is in database, if it is not it adds it
        '''
        rxn_t = tuple([reaction_ID])
        self.all_rxns[reaction_ID] = name, kegg
        self.all_reaction_compound[reaction_ID] = []
        self.all_reaction_compound_store[reaction_ID] = []
        if (reaction_ID, self.mi, genes) not in self.genelist:
            self.genelist.append((reaction_ID, self.mi, genes))
        if (reaction_ID, self.mi, proteins) not in self.proteinlist:
            self.proteinlist.append((reaction_ID, self.mi, proteins))
        if (reaction_ID, self.mi, revers) not in self.modelreactions:
            self.modelreactions.append((reaction_ID, self.mi, revers))
        if (reaction_ID in self.all_rxnreversibility.keys()
                and self.all_rxnreversibility[reaction_ID] != revers):
            self.all_rxnreversibility[reaction_ID] = 'true'
        else:
            self.all_rxnreversibility[reaction_ID] = revers
        for compound in temp_all_rxn_compound:
            total = rxn_t+compound
            self.all_reaction_compound[reaction_ID].append(total)
            self.all_reaction_compound_store[reaction_ID].append(compound)

    def read_metacyc_file(self, BIOCYC_translator, file_name):
        '''
        Reads and parses metacyc SBML file
        '''
        self.count_addedrxns = 0
        self.BIOCYC_translator = BIOCYC_translator
        self.all_compounds = self.DB.get_all_compounds()
        self.modelreactions = []
        self.modelcompounds = []
        self.genelist = []
        self.proteinlist = []
        self.all_reaction_compound = {}
        self.all_rxns = {}
        self.all_cpds = {}
        self.all_rxnreversibility = {}
        self.added_rxns = {}
        self.all_reaction_compound_store = {}
        self.store_compounds = {}
        self.store_compounds_cpt = {}
        if self.inchidb is True:
            self.CT = pit.CompoundTranslator()

        filename = open(file_name)
        strainer = SoupStrainer('model')
        soup = BeautifulSoup(filename, "lxml", parse_only=strainer)
        model = soup.find('model')
        self.mi = model.get('id')
        species_soup = soup.listofspecies
        allspecies = species_soup.findAll('species')
        count = 0
        print ('STATUS: compiling metacyc compounds')
        for c in tqdm(allspecies):
            count += 1
            attr = c.findAll('p')
            inchi_id = None
            biocycID = None
            kegg_id = None
            for a in attr:
                if a.string.startswith('BIOCYC'):
                    biocyc = a.string.split(': ')
                    biocycID = biocyc[1]
                if a.string.startswith('INCHI'):
                    inchi = a.string.split(': ')
                    inchi_id = inchi[1]
                if a.string.startswith('KEGG'):
                    kegg = a.string.split(': ')
                    kegg_id = kegg[1]                
            self.store_compounds_cpt[c.get('id')] = c.get('compartment')+'0'
            self.compound_translator(c.get('id'), biocycID, inchi_id, kegg_id,
                                     c.get('name'), c.get('compartment')+'0')
        reactions_soup = soup.listofreactions
        allreactions = reactions_soup.findAll('reaction')
        print ('STATUS: compiling metacyc reactions')
        for rxn in tqdm(allreactions):
            attr = rxn.findAll('p')
            biocyc_rxn_ID = None
            gene_4_rxn = None
            protein_4_rxn = None
            kegg_4_rxn = None
            for a in attr:
                if a.string.startswith('BIOCYC'):
                    biocyc = a.string.split(': ')
                    biocyc_rxn_ID = biocyc[1]
                if a.string.startswith('GENE'):
                    genes = a.string.split(': ')
                    gene_4_rxn = genes[1]
                if a.string.startswith('EC Number'):
                    proteins = a.string.split(': ')
                    protein_4_rxn = proteins[1]
                if a.string.startswith('KEGG'):
                    kegg = a.string.split(': ')
                    kegg_4_rxn = kegg[1]
            if (gene_4_rxn or protein_4_rxn):
                if biocyc_rxn_ID in self.BIOCYC_translator['rxn']:
                    for krxn in self.BIOCYC_translator['rxn'][biocyc_rxn_ID]:
                        self.retrieve_rxn_info(rxn, krxn, gene_4_rxn, protein_4_rxn, kegg_4_rxn, biocycID=True)
                else:
                    self.retrieve_rxn_info(rxn, rxn.get('id'), gene_4_rxn, protein_4_rxn, kegg_4_rxn, biocycID=False)
            else:
                verbose_print(self.verbose, 'STATUS {} does not have gene or reaction therefore not adding to database'.format(rxn.get('id')))
        print('STATUS: {} number of new metacyc reactions added'.format(self.count_addedrxns))

    def fill_temp_array(self, cpdID, is_prod, stoic, temp_all_rxn_compound):
        '''
        Adds reaction information from metacyc xml file to temporary reaction
        list which later gets inserted into database
        '''
        if (cpdID, is_prod, stoic) not in temp_all_rxn_compound:
            temp_all_rxn_compound.append((cpdID, is_prod, stoic))
        return temp_all_rxn_compound

    def check_db(self, cpdID, kegg_id):
        '''
        Checks database for table original_db_cpdIDs if table exists and cpd ID exists
        the inchi value for that compound is retrieved
        '''
        INCHI_CHECK = None
        KEGG_CHECK = None
        try:
            QC = self.cnx.execute("""SELECT inchi_id from original_db_cpdIDs where ID=?""", (cpdID,))
            try:
                result = QC.fetchone()[0]
                cpdID = result
                INCHI_CHECK = True
            except TypeError:
                pass
        except sqlite3.OperationalError:
            pass
        if not INCHI_CHECK and kegg_id != 'None':
            try:
                QC = self.cnx.execute("""SELECT ID from compound where kegg_id=?""", (kegg_id,))
                try:
                    result = QC.fetchone()[0]
                    cpdID = result
                    KEGG_CHECK = True
                except TypeError:
                    pass
            except sqlite3.OperationalError:
                pass
        return cpdID 

    def get_compounds_4_rxn(self, species):
        '''get reaction compounds'''
        temp_all_rxn_compound = []
        temp_compartment = []
        for comp in species:
            compound_ID = comp['species']
            temp_compartment.append(self.store_compounds_cpt[compound_ID])
            is_prod = (comp.parent.name == "listofproducts")
            stoic = comp['stoichiometry']
            if compound_ID in self.store_compounds:
                temp_all_rxn_compound = self.fill_temp_array(self.store_compounds[compound_ID],
                                                             is_prod, stoic,
                                                             temp_all_rxn_compound)
            else:
                print ('WARNING: compound {} not in dictionary'.format(compound_ID))
                temp_all_rxn_compound = self.fill_temp_array(compound_ID, is_prod,
                                                             stoic, temp_all_rxn_compound)
        return(temp_compartment, temp_all_rxn_compound)

    def retrieve_rxn_info(self, rxn, rxnID, genes, proteins, kegg_ID, biocycID):
        '''
        Parses reaction information from metacyc xml file
        '''
        name = rxn.get('name')
        if kegg_ID:
            QC = self.cnx.execute("""SELECT ID from reaction where kegg_id=?""", (kegg_ID,))
            results = QC.fetchall()
            if results:
                if (results[0][0], self.mi, rxn['reversible'])  not in self.modelreactions:
                    self.modelreactions.append((results[0][0], self.mi, rxn['reversible']))
                    self.proteinlist.append((results[0][0], self.mi, str(proteins)))
                    self.genelist.append((results[0][0], self.mi, str(genes)))
            else:
                temp_compartment, temp_all_rxn_compound = self.get_compounds_4_rxn(rxn.findAll("speciesreference"))
                rxnID = self.retrieve_compartment_4_compartment(temp_compartment, kegg_ID)
                self.multiple_copies_of_rxns(rxnID, rxn, name, genes, proteins, str(kegg_ID), temp_all_rxn_compound)

        elif biocycID:
            temp_compartment, temp_all_rxn_compound = self.get_compounds_4_rxn(rxn.findAll("speciesreference"))
            rxnID = self.retrieve_compartment_4_compartment(temp_compartment, rxnID)
            self.multiple_copies_of_rxns(rxnID, rxn, name, genes, proteins, str(kegg_ID), temp_all_rxn_compound)

        else:
            temp_compartment, temp_all_rxn_compound = self.get_compounds_4_rxn(rxn.findAll("speciesreference"))
            rxnID = self.retrieve_compartment_4_compartment(temp_compartment, rxnID)
            self.rxn_translator(rxnID, temp_all_rxn_compound, rxn['reversible'], name,
                                str(genes), str(proteins), str(kegg_ID))            

    def multiple_copies_of_rxns(self, rxnID, rxn, name, genes, proteins, kegg, temp_all_rxn_compound):
        '''Deals with rxns that have same catalytic enzyme but different substrates'''
        if rxnID not in self.all_reaction_compound:
            self.count_addedrxns += 1
            self.added_rxns[rxnID] = []
            self.added_rxns[rxnID].append(rxnID)
            self.rxn_translator(rxnID, temp_all_rxn_compound, rxn['reversible'], name,
                                str(genes), str(proteins), str(kegg))
        else:
            verbose_print(self.verbose, 'STATUS: Already reaction present {}'.format(rxnID))
            if self.check_reaction_difference(rxnID, temp_all_rxn_compound):
                self.count_addedrxns += 1
                versionrxn = '_v'+str(len(self.added_rxns[rxnID]))
                count_numbofversions = 0
                for kversion in self.added_rxns[rxnID]:
                    if (sorted(self.all_reaction_compound_store[kversion]) ==
                            sorted(temp_all_rxn_compound)):
                        count_numbofversions += 1
                if count_numbofversions == 0:
                    verbose_print(self.verbose, "STATUS: The version of {} reaction added".format(rxnID+versionrxn))
                    self.added_rxns[rxnID].append(rxnID+versionrxn)
                    self.rxn_translator(rxnID+versionrxn, temp_all_rxn_compound, rxn['reversible'],
                                        name, str(genes), str(proteins), str(kegg))        

    def check_reaction_difference(self, rxnID, rxn_compounds):
        '''Checks to see if promiscus mets are only difference between reactions'''
        cpd_ids = [i[0] for i in rxn_compounds]
        cpd_ids = set(cpd_ids)
        
        current_cpd_ids = self.all_reaction_compound_store[rxnID]
        current_cpd_ids = [i[0] for i in current_cpd_ids]
        current_cpd_ids = set(current_cpd_ids)
        diff_set = current_cpd_ids - cpd_ids
        result = diff_set.difference(self.prom_cpds)
        if not result:
            verbose_print(self.verbose, 'STATUS: Reaction only difference may be promiscuous metabolites so not adding {}'.format(rxnID))
            return (False)
        else:
            verbose_print(self.verbose, 'STATUS: Reaction {} difference may include other metabolites than promiscuous metabolites so adding'.format(rxnID))
            return (True)
    def retrieve_compartment_4_compartment(self, temp_compartment, rxnID):
        if len(set(temp_compartment)) == 1:
            rxnID = rxnID+'_'+temp_compartment[0]
        else:
            rxnID = rxnID+'_t0'
        return rxnID