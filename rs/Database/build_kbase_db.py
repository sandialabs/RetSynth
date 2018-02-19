from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'builds tables for sqlite database for kbase file'
import os
import sqlite3
import glob
import re
import urllib2
import httplib
import pubchempy
from copy import deepcopy
from bs4 import BeautifulSoup, SoupStrainer
from Pubchem import pubchem_inchi_translator as pit
from tqdm import tqdm
KEGG = 'http://rest.kegg.jp/'

def parse_data_sbmlfile(inchi, CPD2KEGG, RXN2KEGG, file_name, inchi_pubchem):
    '''
    Open metabolic network file (xml file) and parse information
    '''
    with open(file_name) as f:
        print (file_name+' processing ...')
        strainer = SoupStrainer('model')
        soup = BeautifulSoup(f, "lxml", parse_only=strainer)
        m = soup.find('model')
        mi = m['id']

    modelreactions, rxn_info, all_rxn_cpds, genelist, proteinlist, keggdict = process_reactions(
       soup.listofreactions,  RXN2KEGG, CPD2KEGG, mi, inchi, inchi_pubchem)
    modelcompounds, modelcompounds_allinfo = process_compounds(soup.listofspecies, CPD2KEGG,
                                                               mi, inchi, inchi_pubchem)
    modelcompartments = process_compartments(soup.listofcompartments)
    return(modelcompartments, modelcompounds, modelcompounds_allinfo, modelreactions,
           rxn_info, all_rxn_cpds, genelist, proteinlist, keggdict, mi)

def process_reactions(reaction_soup, RXN2KEGG, CPD2KEGG, mi, inchi, inchi_pubchem):
    '''
    Parse reaction information from metabolic network file (xml file)
    '''
    keggdict = {}
    proteinlist = []
    genelist = []
    modelreactions = []
    rxn_info = {}
    all_rxn_cpds = {}
    model_rxns = reaction_soup.findAll("reaction")
    for rxn in model_rxns:
        if rxn['id'] == 'biomass0':
            rxn['id'] = rxn['id']+'_'+mi
        new_rxn_ID, KEGG_ID = get_KEGG_IDs(rxn['id'], RXN2KEGG)
        rxn_info[new_rxn_ID] = {}
        rxn_info[new_rxn_ID]['name'] = rxn['name']
        rxn_reversibility = rxn.get('reversible', 'true')
        rxn_info[new_rxn_ID]['reversible'] = rxn_reversibility
        rxn_info[new_rxn_ID]['kegg'] = str(KEGG_ID)
        all_rxn_cpds[new_rxn_ID] = []

        for cpd in rxn.findAll("speciesreference"):
            if cpd['species'].endswith('_b'):
                pass
            else:
                new_cpdID, original_KEGG_cpdID = get_KEGG_IDs(cpd['species'], CPD2KEGG)
                if inchi:
                    cpdID = inchi_pubchem.get(new_cpdID, new_cpdID)
                else:
                    cpdID = new_cpdID
                is_prod = (cpd.parent.name == "listofproducts")
                stoic = cpd.get('stoichiometry', 1)

                if (new_rxn_ID, cpdID, is_prod, stoic) not in all_rxn_cpds[new_rxn_ID]:
                    all_rxn_cpds[new_rxn_ID].append((new_rxn_ID, cpdID, is_prod, stoic))
        '''Retrieve gene associations'''
        try:
            associations = [a.get_text().strip() for a in rxn.notes.findChildren()]
            for a in associations:
                a = str(a)
                typeassociation, value = a.split(':')
                if typeassociation == 'GENE_ASSOCIATION':
                    if value != '' and value != 'Unknown':
                        genelist.append((new_rxn_ID, mi, value))
                if typeassociation == 'PROTEIN_CLASS':
                    if value != '' and value != 'Unknown':
                        proteinlist.append((new_rxn_ID, mi, value))

        except AttributeError:
            try:
                associations = rxn.p.string
                associations = re.sub('\n', '', associations)
                associations = re.sub('\s+', '', associations)
                typeassociation, value = associations.split(':')
                if typeassociation == 'GENE_ASSOCIATION':
                    if value != '' and value != 'Unknown':
                        genelist.append((new_rxn_ID, mi, value))
            except AttributeError:
                pass
        ''' Get reaction id and reversibility'''
        modelreactions.append((new_rxn_ID, mi, rxn.get('reversible', 'true')))

    return(modelreactions, rxn_info, all_rxn_cpds, genelist, proteinlist, keggdict)

def process_compounds(speciessoup, CPD2KEGG, mi, inchi, inchi_pubchem):
    '''
    Parse compound information from metabolic network file (xml file)
    '''
    modelcompounds = []
    modelcompounds_allinfo = []
    model_cpds = speciessoup.findAll("species")
    for m in model_cpds:
        new_cpdID, original_KEGG_cpdID = get_KEGG_IDs(m['id'], CPD2KEGG)

        try:
            if m['boundarycondition'] == 'false':
                if inchi:
                    modelcompounds.append((inchi_pubchem.get(new_cpdID, new_cpdID), mi))
                    modelcompounds_allinfo.append((inchi_pubchem.get(new_cpdID, new_cpdID), m['name'], m['compartment'], str(original_KEGG_cpdID)))
                else:
                    modelcompounds.append((new_cpdID, mi))
                    modelcompounds_allinfo.append((new_cpdID, m['name'], m['compartment'], str(original_KEGG_cpdID)))
        except KeyError:
            if inchi:
                modelcompounds.append((inchi_pubchem.get(new_cpdID, new_cpdID), mi))
                modelcompounds_allinfo.append((inchi_pubchem.get(new_cpdID, new_cpdID),
                                               m['name'], m['compartment'], str(original_KEGG_cpdID)))
            else:
                modelcompounds.append((new_cpdID, mi))
                modelcompounds_allinfo.append((new_cpdID, m['name'], m['compartment'], str(original_KEGG_cpdID)))
    return(modelcompounds, modelcompounds_allinfo)

def process_compartments(compartmentsoup):
    '''Get compartments in xml file'''
    compartments_all = []
    compartments = compartmentsoup.findAll("compartment")
    for compartment in compartments:
        if (compartment['id'], compartment['name']) not in compartments_all:
            compartments_all.append((compartment['id'], compartment['name']))
    return compartments_all

def get_KEGG_IDs(ID, KEGGdict):
    '''Retrieve KEGG IDs'''
    compartment = get_compartment_info(ID)
    ID_copy = deepcopy(ID)
    ID_kbase = re.sub(compartment+'$', '', ID_copy)
    try:
        KEGG_ID = KEGGdict[ID_kbase]
        original_KEGG_ID = str(deepcopy(KEGG_ID))
        if not KEGG_ID:
            KEGG_ID = ID
        else:
            KEGG_ID = KEGG_ID+compartment
    except KeyError:
        KEGG_ID = ID
        original_KEGG_ID = 'None'
    return(KEGG_ID, original_KEGG_ID)

def open_translation_file(file_name):
    '''opens and stores KEGG translation files '''
    dictionary = {}
    with open(file_name) as fin:
        for line in fin:
            line = line.strip()
            larray = line.split('\t')
            try:
                KEGGIDS = larray[1].split('|')
                dictionary[larray[0]] = KEGGIDS[0]
            except IndexError:
                dictionary[larray[0]] = None
    return dictionary 

def BuildKbase(sbml_dir, kbase2keggCPD_translate_file, kbase2keggRXN_translate_file, inchi, DBpath, rxntype='bio'):
    '''
    Inserts values from metabolic networks, xml into sqlite database
    '''


    '''Load KEGG cpd and rxn conversions'''
    CPD2KEGG = open_translation_file(kbase2keggCPD_translate_file)
    RXN2KEGG = open_translation_file(kbase2keggRXN_translate_file)

    inchi_pubchem = {}
    total_set = set()
    cnx = sqlite3.connect(DBpath)
    cnx.execute("PRAGMA synchronous = OFF")
    cnx.execute("PRAGMA journal_mode = OFF")
    cnx.commit()
    sbml_files = glob.glob(os.path.join(sbml_dir, '*'))
    sbml_files_individual = [os.path.basename(file)
                             for file in glob.glob(os.path.join(sbml_dir, '*'))]

    if inchi:
        print ('Retrieving inchi values for compounds ...')
        CT = pit.CompoundTranslator()
        for filename in sbml_files:
            inchi_pubchem, total_set = get_inchi_values(filename, inchi_pubchem, total_set, CPD2KEGG, CT)

    args = [(i, CPD2KEGG, RXN2KEGG, DBpath, inchi, inchi_pubchem, c, sbml_files_individual[c], rxntype)
            for c, i in enumerate(sbml_files)]

    print ('STATUS: Loading file information into database ...')
    for arg in args:
        load_file_info_2_db(arg)
    print ('STATUS: Finished loading file information into database')
    print ('STATUS: Identifying metabolic clusters...')
    retrieve_metabolic_clusters(DBpath)
    print ('STATUS: Finished identifying metabolic clusters...')
    print ('STATUS: Removing duplicate entries ...')
    if inchi:
        cnx.execute("""DELETE FROM original_db_cpdIDs WHERE rowid NOT IN\
                    (SELECT min(rowid) from original_db_cpdIDs group by ID)""")

    cnx.execute("""DELETE FROM compound WHERE rowid NOT IN\
                (SELECT min(rowid) from compound group by ID)""")

    cnx.execute("""DELETE FROM reaction WHERE rowid NOT IN\
                (SELECT min(rowid) from reaction group by ID)""")

    cnx.execute("""DELETE FROM compartments WHERE rowid NOT IN\
                (SELECT min(rowid) from compartments group by ID)""")

    cnx.execute("""DELETE FROM reaction_reversibility WHERE rowid  NOT IN\
                (SELECT min(rowid) from reaction_reversibility\
                 group by reaction_ID, is_reversible)""")

    cnx.execute("""DELETE FROM reaction_compound WHERE rowid NOT IN\
                 (SELECT min(rowid) from reaction_compound group by\
                  reaction_ID,cpd_ID,is_prod)""")

    Q = cnx.execute("""SELECT reaction_ID FROM reaction_reversibility""")
    results = Q.fetchall()
    for result in results:
        Q = cnx.execute("""SELECT reaction_ID FROM reaction_reversibility WHERE reaction_ID = ?""",
                        (result[0],))
        results1 = Q.fetchall()
        if len(results1) == 2:
            cnx.execute("""DELETE FROM reaction_reversibility\
             WHERE reaction_ID = ? AND is_reversible = ?""", (result[0], 'false'))
        Q = cnx.execute("""SELECT filenum from reaction_compound where reaction_ID = ?""",
                        (result[0],))
        hits = list(set(Q.fetchall()))
        if len(hits) == 2:
            cnx.execute("""DELETE FROM reaction_compound WHERE reaction_ID = ? and filenum = ?""",
                        (result[0], hits[0][0]))

            Q = cnx.execute("""SELECT is_reversible FROM reaction_reversibility \
            	             WHERE reaction_ID=?""", (result[0],))

            hits = Q.fetchone()
            if hits[0] == 'false':
                command = """UPDATE reaction_reversibility SET\
                			 is_reversible = 'true' WHERE reaction_ID = ?"""
                cnx.execute(command, (result[0],))
    cnx.commit()
    print ('STATUS: Finished removing duplicate entries')

def extract_KEGG_data(url):
    '''Extract Kegg db info'''
    try:
        data = urllib2.urlopen(url).read()
        darray = data.split('\n')
        return(darray)
    except urllib2.HTTPError:
        return(None)

def kegg2pubcheminchi(cpd):
    '''Using KEGG ID  to get inchi '''
    darray = extract_KEGG_data(KEGG+'get/'+cpd)
    if darray:
        for value in darray:
            array = value.split()
            if 'PubChem:' in array:
                index = array.index('PubChem:')
                sid = array[index+1]
                substance = pubchempy.Substance.from_sid(sid)
                try:
                    substance_cids = substance.cids
                    if substance_cids:
                        try:
                            compounds = pubchempy.get_compounds(substance_cids[0])
                            if compounds:
                                return compounds[0].inchi
                            else:
                                return None
                        except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                            return None
                except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                    return None                
            else:
                return None
    else:
        return None

def retrieve_exact_inchi_values(m, total_set, inchi_pubchem, CPD2KEGG, CT):
    '''Retrieve InChI values'''
    if m['id'] not in total_set:
        new_KEGG_ID, original_KEGG_ID = get_KEGG_IDs(m['id'], CPD2KEGG)
        compart_info = get_compartment_info(m['id'])
        if original_KEGG_ID != 'None':
            inchi = kegg2pubcheminchi(original_KEGG_ID)
            if not inchi:
                inchi, iupac_name = CT.translate(m['name'])
        else:
            inchi, iupac_name = CT.translate(m['name'])

        if inchi is not None:
            inchi_pubchem[new_KEGG_ID] = inchi+compart_info
            total_set.add(m['id'])
        else:
            total_set.add(m['id'])
    return (inchi_pubchem, total_set)

def get_inchi_values(file_name, inchi_pubchem, total_set, CPD2KEGG, CT):
    '''
    Retrieve InChI values for compounds in metabolic networks
    '''
    with open(file_name, 'r') as f:
        print ('Retrieving  InChi translation for compounds in '+file_name)
        strainer = SoupStrainer('model')
        soup = BeautifulSoup(f, "lxml", parse_only=strainer)
        speciessoup = soup.listofspecies
        model_cpds = speciessoup.findAll("species")
        for m in tqdm(model_cpds):
            try:
                if m['boundarycondition'] == 'false':
                    inchi_pubchem, total_set = retrieve_exact_inchi_values(m, total_set, inchi_pubchem, CPD2KEGG, CT)

            except KeyError:
                inchi_pubchem, total_set = retrieve_exact_inchi_values(m, total_set, inchi_pubchem, CPD2KEGG, CT)
    
    return(inchi_pubchem, total_set)

def get_compartment_info(ID):
    '''
    Retrieve compartment information for compound
    '''

    compart_info = ''
    match = re.search('_\w{1}\d*$', ID)
    if match is not None:
        compart_info = match.group(0)
    else:
        pass
    return compart_info

def load_file_info_2_db(args):
    '''
    Open and insert information from metabolic networks
    (xml files) into database
    '''
    filename = args[0]
    CPD2KEGG = args[1]
    RXN2KEGG = args[2]
    DBpath = args[3]
    inchi = args[4]
    inchi_pubchem = args[5]
    filenum = args[6]
    rawfilename = args[7]
    rxntype = args[8]
    (modelcompartments, modelcompounds, modelcompounds_allinfo, modelreactions, rxn_info,
     all_rxn_cpds, genelist, proteinlist, keggdict, mi) = parse_data_sbmlfile(inchi, CPD2KEGG, RXN2KEGG,
                                                                              filename, inchi_pubchem)
    insert_individual_model_results_2_db(DBpath, modelcompounds, modelreactions,
                                         genelist, proteinlist, mi, rawfilename)
    insert_comprehensive_model_results_2_db(DBpath, modelcompartments, modelcompounds_allinfo, rxn_info,
                                            all_rxn_cpds, keggdict, inchi, inchi_pubchem, filenum, rxntype)

def insert_individual_model_results_2_db(DBpath, modelcompounds, modelreactions,
                                         genelist, proteinlist, mi, filename):
    '''
    Inserts individual model information into tables
    '''
    cnx = sqlite3.connect(DBpath)
    cnx.executemany("INSERT INTO model_compound VALUES (?,?)", modelcompounds)
    cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", modelreactions)

    cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", genelist)
    cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", proteinlist)

    cnx.execute("INSERT INTO model VALUES (?,?)", (mi, filename))
    cnx.commit()

def insert_comprehensive_model_results_2_db(DBpath, modelcompartments, modelcompounds_allinfo, rxn_info,
                                            all_rxn_cpds, keggdict, inchi, inchi_pubchem, filenum, rxntype):
    '''
    Inserts comprehensive compound and reaction information into tables
    '''
    cnx = sqlite3.connect(DBpath)
    cnx.executemany("INSERT INTO compartments VALUES (?,?)", modelcompartments)

    def add_rxns_2_db(cnx, rxn_info, all_rxn_cpds):
        '''add reaciton to database'''
        for rxn in rxn_info:
            Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", (rxn,))
            result = Q.fetchone()
            if result is None:
                cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)",
                            (rxn, rxn_info[rxn]['reversible']))
                cnx.execute("INSERT INTO reaction VALUES (?,?,?, ?)", (rxn, rxn_info[rxn]['name'], rxn_info[rxn]['kegg'], rxntype))
                cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)", (rxn, rxn_info[rxn]['reversible']))
                for rxn_compound in all_rxn_cpds[rxn]:
                    cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                     (rxn_compound[0], rxn_compound[1],
                                      rxn_compound[2], rxn_compound[3], filenum))
            else:
                Q = cnx.execute("""SELECT is_reversible FROM reaction_reversibility \
                                     WHERE reaction_ID = ?""", (rxn,))
                result = Q.fetchone()
                if result[0] != rxn_info[rxn]['reversible'] and result[0] == 'false':
                    cnx.execute("""UPDATE reaction_reversibility SET is_reversible = 'true' WHERE reaction_ID = ?""",(rxn,))
        cnx.commit()

    def add_cpds_2_db(cnx, modelcompounds_allinfo):
        '''add compounds to database'''
        for cpd in modelcompounds_allinfo:
            Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?", (cpd[0],))
            result = Q.fetchone()
            if result is None:
                cnx.execute("""INSERT INTO compound VALUES (?,?,?,?)""",
                                 (cpd[0], cpd[1], cpd[2], cpd[3]))
        cnx.commit()
    add_rxns_2_db(cnx, rxn_info, all_rxn_cpds)
    add_cpds_2_db(cnx, modelcompounds_allinfo)

    if inchi is True:
        for ID, inchi_id in inchi_pubchem.iteritems():
            cnx.execute("""INSERT INTO original_db_cpdIDs VALUES (?,?)""", (ID, inchi_id))
    cnx.commit()

def retrieve_metabolic_clusters(DBpath):
    '''
    Identifies and inserts metabolic clusters
    (organisms with the same compounds and reactions) into database
    '''
    cnx = sqlite3.connect(DBpath)
    metabolic_clusters = {}
    cluster_org = {}
    Q = cnx.execute("SELECT ID FROM model")
    results = Q.fetchall()
    for result in results:
        Q = cnx.execute("SELECT cpd_ID FROM model_compound WHERE model_ID = ?",
                        (result[0],))
        resultscpds = Q.fetchall()
        org_cpds = [cpd[0] for cpd in resultscpds]
        Q = cnx.execute("SELECT reaction_ID FROM model_reaction WHERE model_ID = ?", (result[0],))
        resultsrxns = Q.fetchall()
        org_rxns = []
        for rxn in resultsrxns:
            if rxn[0].startswith('biomass'):
                pass
            else:
                org_rxns.append(rxn[0])
        metabolic_clusters, cluster_org = get_metabolic_clusters(sorted(org_cpds),
                                                                 sorted(org_rxns), result[0],
                                                                 metabolic_clusters, cluster_org)

    Q = cnx.execute("DELETE FROM cluster")
    for key, value in cluster_org.iteritems():
        for v in value:
            Q = cnx.execute("SELECT cluster_num FROM cluster WHERE cluster_num= ? AND ID=?",
                            (key, v))
            result = Q.fetchone()
            if result is None:
                cnx.execute("INSERT INTO cluster VALUES (?,?)", (key, v))
    cnx.commit()

def get_metabolic_clusters(org_cpds, org_rxns, model_id,
                           metabolic_clusters, cluster_org):
    '''
    Retrieves metabolic cluster info
    '''
    if len(metabolic_clusters) > 0:
        count_keys = 0
        for key in metabolic_clusters:
            if (org_cpds == metabolic_clusters[key]['comps'] and
                    org_rxns == metabolic_clusters[key]['rxns']):
                cluster_org[key].append(model_id)
                break
            else:
                count_keys += 1
        if len(metabolic_clusters) == count_keys:
            num_keys = len(metabolic_clusters)
            metabolic_clusters[num_keys+1] = {}
            metabolic_clusters[num_keys+1]['comps'] = org_cpds
            metabolic_clusters[num_keys+1]['rxns'] = org_rxns
            cluster_org[num_keys+1] = []
            cluster_org[num_keys+1].append(model_id)
    else:
        metabolic_clusters[1] = {}
        metabolic_clusters[1]['rxns'] = org_rxns
        metabolic_clusters[1]['comps'] = org_cpds
        cluster_org[1] = []
        cluster_org[1].append(model_id)
    return(metabolic_clusters, cluster_org)
