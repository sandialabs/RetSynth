from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'builds tables for sqlite database'
import os
import sqlite3
import glob
import re
from bs4 import BeautifulSoup, SoupStrainer
from Pubchem import pubchem_inchi_translator as pit
from tqdm import tqdm
def parse_data_sbmlfile(inchi, file_name, inchi_pubchem):
    '''
    Open metabolic network file (xml file) and parse information
    '''
    with open(file_name) as f:
        print (file_name+' processing ...')
        strainer = SoupStrainer('model')
        soup = BeautifulSoup(f, "lxml", parse_only=strainer)
        m = soup.find('model')
        mi = m['id']

    modelreactions, rxn_info, all_rxn_cpds, genelist, proteinlist = process_reactions(
        soup.listofreactions, mi, inchi, inchi_pubchem)
    modelcompounds, modelcompounds_allinfo = process_compounds(soup.listofspecies,
                                                               mi, inchi, inchi_pubchem)
    modelcompartments = process_compartments(soup.listofcompartments)
    return(modelcompartments, modelcompounds, modelcompounds_allinfo, modelreactions,
           rxn_info, all_rxn_cpds, genelist, proteinlist, mi)

def process_reactions(reaction_soup, mi, inchi, inchi_pubchem):
    '''
    Parse reaction information from metabolic network file (xml file)
    '''
    proteinlist = []
    genelist = []
    modelreactions = []
    rxn_info = {}
    all_rxn_cpds = {}
    model_rxns = reaction_soup.findAll("reaction")
    for rxn in model_rxns:
        if rxn['id'] == 'biomass0':
            rxn['id'] = rxn['id']+'_'+mi
        rxn_info[rxn['id']] = {}
        rxn_info[rxn['id']]['name'] = rxn['name']
        rxn_reversibility = rxn.get('reversible', 'true')
        rxn_info[rxn['id']]['reversible'] = rxn_reversibility
        all_rxn_cpds[rxn['id']] = []

        for cpd in rxn.findAll("speciesreference"):
            if cpd['species'].endswith('_b'):
                pass
            else:
                if inchi:
                    cpdID = inchi_pubchem.get(cpd['species'], cpd['species'])
                else:
                    cpdID = cpd['species']
                is_prod = (cpd.parent.name == "listofproducts")
                stoic = cpd.get('stoichiometry', 1)

                if (rxn['id'], cpdID, is_prod, stoic) not in all_rxn_cpds[rxn['id']]:
                    all_rxn_cpds[rxn['id']].append((rxn['id'], cpdID, is_prod, stoic))
        '''Retrieve gene associations'''
        try:
            associations = [a.get_text().strip() for a in rxn.notes.findChildren()]
            for a in associations:
                a = str(a)
                typeassociation, value = a.split(':')
                if typeassociation == 'GENE_ASSOCIATION':
                    if value != '' and value != 'Unknown':
                        genelist.append((rxn['id'], mi, value))
                if typeassociation == 'PROTEIN_CLASS':
                    if value != '' and value != 'Unknown':
                        proteinlist.append((rxn['id'], mi, value))
        except AttributeError:
            try:
                associations = rxn.p.string
                associations = re.sub('\n', '', associations)
                associations = re.sub('\s+', '', associations)
                typeassociation, value = associations.split(':')
                if typeassociation == 'GENE_ASSOCIATION':
                    if value != '' and value != 'Unknown':
                        genelist.append((rxn['id'], mi, value))
            except AttributeError:
                pass
        ''' Get reaction id and reversibility'''
        modelreactions.append((rxn['id'], mi, rxn.get('reversible', 'true')))

    return(modelreactions, rxn_info, all_rxn_cpds, genelist, proteinlist)

def process_compounds(speciessoup, mi, inchi, inchi_pubchem):
    '''
    Parse compound information from metabolic network file (xml file)
    '''
    modelcompounds = []
    modelcompounds_allinfo = []
    model_cpds = speciessoup.findAll("species")
    for m in model_cpds:
        try:
            if m['boundarycondition'] == 'false':
                if inchi:
                    modelcompounds.append((inchi_pubchem.get(m['id'], m['id']), mi))
                    modelcompounds_allinfo.append((inchi_pubchem.get(m['id'], m['id']),
                                                   m['name'], m['compartment']))
                else:
                    modelcompounds.append((m['id'], mi))
                    modelcompounds_allinfo.append((m['id'], m['name'], m['compartment']))
        except KeyError:
            if inchi:
                modelcompounds.append((inchi_pubchem.get(m['id'], m['id']), mi))
                modelcompounds_allinfo.append((inchi_pubchem.get(m['id'], m['id']),
                                               m['name'], m['compartment']))
            else:
                modelcompounds.append((m['id'], mi))
                modelcompounds_allinfo.append(m['id'], m['name'], m['compartment'])
    return(modelcompounds, modelcompounds_allinfo)

def process_compartments(compartmentsoup):
    '''Get compartments in xml file'''
    compartments_all = []
    compartments = compartmentsoup.findAll("compartment")
    for compartment in compartments:
        if (compartment['id'], compartment['name']) not in compartments_all:
            compartments_all.append((compartment['id'], compartment['name']))
    return compartments_all

def BuildTables(sbml_dir, inchi, DBpath, rxntype='bio'):
    '''
    Inserts values from metabolic networks, xml into sqlite database
    '''
    inchi_pubchem = {}
    no_inchi = set()
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
            inchi_pubchem, no_inchi = get_inchi_values(filename, inchi_pubchem, no_inchi, CT)

    args = [(i, DBpath, inchi, inchi_pubchem, c, sbml_files_individual[c], rxntype)
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

def get_inchi_values(file_name, inchi_pubchem, no_inchi, CT):
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
            total_set = set(inchi_pubchem.keys())
            total_set.update(set(no_inchi))
            try:
                if m['boundarycondition'] == 'false':
                    if m['id'] not in total_set:
                        compart_info = get_compartment_info(m['id'])
                        inchi, iupac_name = CT.translate(m['name'])
                        if inchi is not None:
                            inchi_pubchem[m['id']] = inchi+compart_info
                        else:
                            no_inchi.add(m['id'])

            except KeyError:
                if m['id'] not in total_set:
                    compart_info = get_compartment_info(m['id'])
                    inchi, iupac_name = CT.translate(m['name'])
                    if inchi is not None:
                        inchi_pubchem[m['id']] = inchi+compart_info
                    else:
                        no_inchi.add(m['id'])

    return(inchi_pubchem, no_inchi)

def get_compartment_info(compoundID):
    '''
    Retrieve compartment information for compound
    '''

    compart_info = ''
    match = re.search('_\w{1}\d*$', compoundID)
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
    DBpath = args[1]
    inchi = args[2]
    inchi_pubchem = args[3]
    filenum = args[4]
    rawfilename = args[5]
    rxntype = args[6]
    (modelcompartments, modelcompounds, modelcompounds_allinfo, modelreactions, rxn_info,
     all_rxn_cpds, genelist, proteinlist, mi) = parse_data_sbmlfile(inchi,
                                                                    filename, inchi_pubchem)
    insert_individual_model_results_2_db(DBpath, modelcompounds, modelreactions,
                                         genelist, proteinlist, mi, rawfilename)
    insert_comprehensive_model_results_2_db(DBpath, modelcompartments, modelcompounds_allinfo, rxn_info,
                                            all_rxn_cpds, inchi, inchi_pubchem, filenum, rxntype)

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
                                            all_rxn_cpds, inchi, inchi_pubchem, filenum, rxntype):
    '''
    Inserts comprehensive compound and reaction information into tables
    '''
    cnx = sqlite3.connect(DBpath)
    cnx.executemany("INSERT INTO compartments VALUES (?,?)", modelcompartments)
    cnx.executemany("INSERT INTO compound VALUES (?,?,?)", modelcompounds_allinfo)
    for rxn in rxn_info:
        '''Name info'''
        cnx.execute("INSERT INTO reaction VALUES (?,?,?)", (rxn, rxn_info[rxn]['name'], rxntype))
        cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)",
                    (rxn, rxn_info[rxn]['reversible']))
    for key, rxn_cpds in all_rxn_cpds.iteritems():
        for rxn_cpd in rxn_cpds:
            cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                        (rxn_cpd[0], rxn_cpd[1], rxn_cpd[2], rxn_cpd[3], filenum))

    if inchi is True:
        inchi_id = inchi_pubchem.items()
        cnx.executemany("INSERT INTO original_db_cpdIDs VALUES (?,?)", inchi_id)
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
