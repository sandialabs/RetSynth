from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Build metabolic db using ATLAS db'
KEGG = 'http://rest.kegg.jp/'

from multiprocessing import Process, Queue
import re
import glob
import os
import sqlite3
import urllib2
import httplib
import pubchempy
from tqdm import tqdm

def build_atlas(atlas_dir, DBPath, inchidb, processors, rxntype='bio'):
    '''Add atlas database to RSA metabolic database'''
    atlas_files = glob.glob(os.path.join(atlas_dir, '*'))

    rxn_keggids, rxn_atlas = open_atlas_files(atlas_files)
    (cnx, reactions, reaction_reversible,
     model_reaction, reaction_protein,
     reaction_genes, model_compound,
     compound, reaction_compound, original_db_cpd_new) = fill_arrays_4_db(rxn_keggids, rxn_atlas,
                                                                          DBPath, inchidb,
                                                                          processors, rxntype)
    fill_database(cnx, reactions, reaction_reversible, model_reaction,
                  reaction_protein, reaction_genes, model_compound,
                  compound, reaction_compound, original_db_cpd_new)

def fill_dictionary(larray, dictionary, KEGGID=False):
    '''Fill dictionaries with ATLAS reactions'''
    if KEGGID is True:
        ID = larray[1]
    else:
        ID = larray[0]
    dictionary[ID] = {}
    larray[2] = re.sub('\s+', '', larray[2])
    dictionary[ID]['reactionformula'] = larray[2]
    rxnrules = larray[3].split('|')
    dictionary[ID]['reactionrule'] = rxnrules
    if larray[6] != '' and larray[6] != 'TRUE' and larray[6] != 'FALSE':
        most_similar = larray[6].split('/')
        try:
            ECenzymes = most_similar[1].split('|')
            if ECenzymes:
                dictionary[ID]['mostsim'] = ECenzymes
        except IndexError:
            pass
    return dictionary

def fill_dictionary_atlasbiochem(larray, dictionary):
    '''Fill dictionaries with ATLAS reactions'''
    ID = larray[0]
    if ID not in dictionary:
        dictionary[ID] = {}
        larray[1] = re.sub('\s+', '', larray[1])
        dictionary[ID]['reactionformula'] = larray[1]
        if larray[2] != '':
            dictionary[ID]['name'] = larray[2]
        rxnrules = larray[4].split(' | ')
        dictionary[ID]['reactionrule'] = rxnrules
        if larray[3] != '':
            ECenzymes = larray[3].split(' | ')
            if ECenzymes:
                dictionary[ID]['mostsim'] = ECenzymes
    else:
        if 'name' not in dictionary[ID]:
            if larray[2] != '':
                dictionary[ID]['name'] = larray[2]
    return dictionary

def open_atlas_files(atlas_files):
    '''open ATLAS files'''
    print ('STATUS: opening and retrieving information from ATLAS files')
    rxn_keggids = {}
    rxn_atlas = {}
    for filename in atlas_files:
        if filename.endswith('ATLAS-FULL.csv'):
            print ('STATUS: opening {}'.format(filename))
            with open(filename) as fin:
                header = fin.readline()
                for line in fin:
                    line = line.strip()
                    larray = line.split(',')
                    if larray[1] != '':
                        rxn_keggids = fill_dictionary(larray, rxn_keggids, KEGGID=True)
                    else:
                        rxn_atlas = fill_dictionary(larray, rxn_atlas)
        else:
            print ('STATUS: opening {}'.format(filename))
            with open(filename) as fin:
                header = fin.readline()
                for line in fin:
                    line = line.strip()
                    line = re.sub('"', '', line)
                    larray = line.split(',')
                    rxn_keggids = fill_dictionary_atlasbiochem(larray, rxn_keggids)
    return(rxn_keggids, rxn_atlas)

def fill_arrays_4_db(rxn_keggids, rxn_atlas, DBPath, inchidb, processors, rxntype):
    '''fill arrays with ATLAS reactions to add to database'''
    cnx = sqlite3.connect(DBPath)
    reaction_compound = []
    compound = []
    reactions = []
    reaction_reversible = []
    reaction_protein = []
    reaction_genes = []
    model_compound = []
    model_reaction = []
    original_db_cpd_new = []
    original_db_cpd_current = []
    if inchidb is True:
        Q = cnx.execute("""SELECT * from original_db_cpdIDs""")
        hits = Q.fetchall()
        original_db_cpd_current = list(set(hits))

    Q = cnx.execute("""SELECT kegg_id from reaction""")
    hits = Q.fetchall()
    dbrxns = [i[0] for i in hits]

    Q = cnx.execute("""SELECT ID from reaction""")
    hits = Q.fetchall()
    dbrxns_id = [i[0] for i in hits]
    dbrxns = list(set(dbrxns+dbrxns_id))

    Q = cnx.execute("""SELECT kegg_id from compound""")
    hits = Q.fetchall()
    dbcpds = [i[0] for i in hits]

    Q = cnx.execute("""SELECT ID from compound""")
    hits = Q.fetchall()
    dbcpds_id = [i[0] for i in hits]
    dbcpds = list(set(dbcpds+dbcpds_id))

    currentcpds = {}
    A = cnx.execute('SELECT * from compound')
    results = A.fetchall()
    for i in results:
        currentcpds.setdefault(i[3], {})
        currentcpds[i[3]].setdefault(i[2], i[0])
    
    rxns = []
    print ('STATUS: Retrieving ATLAS rxns with kegg ids')
    for rxn in tqdm(rxn_keggids):
        if rxn not in dbrxns or rxn+'_c0' not in dbrxns:
            try:
                name = rxn_keggids[rxn]['name']
            except KeyError:
                name = 'None'
            try:
                proteins = rxn_keggids[rxn]['mostsim']
            except KeyError:
                proteins = rxn_keggids[rxn]['reactionrule']
            try:
                rxnformula = rxn_keggids[rxn]['reactionformula']
            except KeyError:
                print ('WARNING: no formula for rxn {}'.format(rxn))
                rxnformula = ''
            rxns.append((rxn, name, rxnformula, proteins, True))

    print ('STATUS: Retrieving ATLAS rxns without kegg id')
    for rxn in tqdm(rxn_atlas):
        if rxn not in dbrxns or rxn+'_c0' not in dbrxns:
            try:
                proteins = rxn_atlas[rxn]['mostsim']
            except KeyError:
                proteins = rxn_atlas[rxn]['reactionrule']
            try:
                rxnformula = rxn_atlas[rxn]['reactionformula']
            except KeyError:
                print ('WARNING: no formula for rxn {}'.format(rxn))
                rxnformula = ''
            rxns.append((rxn, 'None', rxnformula, proteins, False))
    output_queue = Queue()
    rxns_processors = [rxns[i:i+processors]
                     for i in range(0, len(rxns), processors)]
    print ('STATUS: Processing all reactions')
    for rxns_processor in tqdm(rxns_processors):
        processes = []
        for rxninfo in rxns_processor:
            processes.append(Process(target=process_reactions,
                                     args=(rxninfo, currentcpds, dbcpds, original_db_cpd_current,
                                           inchidb, rxntype, output_queue)))
        for p in processes:
            p.start()
        for p in processes:
            result_tuples = output_queue.get()
            reactions.extend(result_tuples[0])
            reaction_reversible.extend(result_tuples[1])
            model_reaction.extend(result_tuples[2])
            reaction_protein.extend(result_tuples[3])
            reaction_genes.extend(result_tuples[4])
            model_compound.extend(result_tuples[5])
            compound.extend(result_tuples[6])
            reaction_compound.extend(result_tuples[7])
            original_db_cpd_new.extend(result_tuples[8])   
    reactions = list(set(reactions))
    model_reaction = list(set(model_reaction))
    model_compound = list(set(model_compound))
    compound = list(set(compound))
    reaction_reversible = list(set(reaction_reversible))
    if original_db_cpd_new:
        original_db_cpd_new = list(set(original_db_cpd_new))

    return (cnx, reactions, reaction_reversible, model_reaction, reaction_protein,
            reaction_genes, model_compound, compound, reaction_compound, original_db_cpd_new)

def process_reactions(rxninfo, currentcpds, dbcpds, original_db_cpd_current, inchidb, rxntype, output_queue):
    '''Process ATLAS reactions'''
    reactions_temp = []
    reaction_reversible_temp = []
    model_reaction_temp = []
    reaction_genes_temp = []
    reaction_protein_temp = []
    reaction_compound_temp = []
    model_compound_temp = []
    compound_temp = []
    original_db_cpd_temp = []
    if rxninfo[4] is True:
        reactions_temp.append((rxninfo[0]+'_c0', rxninfo[1], rxninfo[0], rxntype))
    else:
        reactions_temp.append((rxninfo[0]+'_c0', rxninfo[1], 'None', rxntype))
    reaction_reversible_temp.append((rxninfo[0]+'_c0', 'true'))
    model_reaction_temp.append((rxninfo[0]+'_c0', 'ATLAS', 'true'))
    reaction_genes_temp.append((rxninfo[0]+'_c0', 'ATLAS', 'None'))
    
    for protein in rxninfo[3]:
        reaction_protein_temp.append((rxninfo[0]+'_c0', 'ATLAS', protein))

    reactionformula = rxninfo[2]
    reactionarray = reactionformula.split('<=>')
    reactants = reactionarray[0].split('+')
    products = reactionarray[1].split('+')
    if reactants and products:
        for reactant in reactants:
            (model_compound_temp,
             compound_temp,
             reaction_compound_temp,
             original_db_cpd_temp) = process_substrates(rxninfo[0]+'_c0', reactant, 0,
                                                        currentcpds, dbcpds, inchidb,
                                                        original_db_cpd_current,
                                                        model_compound_temp,
                                                        reaction_compound_temp,
                                                        model_reaction_temp,
                                                        compound_temp,
                                                        original_db_cpd_temp)
        for product in products:
            (model_compound_temp,
             compound_temp,
             reaction_compound_temp,
             original_db_cpd_temp) = process_substrates(rxninfo[0]+'_c0', product, 1,
                                                        currentcpds, dbcpds, inchidb,
                                                        original_db_cpd_current,
                                                        model_compound_temp,
                                                        reaction_compound_temp,
                                                        model_reaction_temp,
                                                        compound_temp,
                                                        original_db_cpd_temp)
    else:
        print ('WARNING: {} does not have reactants {} and/or products {}'.format(rxninfo[0], reactants, products))

    output_queue.put((reactions_temp, reaction_reversible_temp, model_reaction_temp,
                      reaction_protein_temp, reaction_genes_temp, model_compound_temp,
                      compound_temp, reaction_compound_temp, original_db_cpd_temp))

def process_substrates(rxn, cpd, is_prod, currentcpds, dbcpds, inchidb, original_db_cpd_current,
                       model_compound_temp, reaction_compound_temp, model_reaction_temp,
                       compound_temp, original_db_cpd_temp):
    '''Process ATLAS compounds'''

    match = re.search('^\(\d+\)', cpd)
    if match is not None:
        stoich = match.group(0)
        stoich = re.sub('\(|\)', '', stoich)
        stoich = int(stoich)
    else:
        stoich = 1
    cpd = re.sub('^\(\d+\)', '', cpd)
    if cpd in currentcpds and 'c0' in currentcpds[cpd]:
        db_cpd_id = currentcpds[cpd]['c0']
        reaction_compound_temp.append((rxn, db_cpd_id, is_prod, stoich, 0))
        model_compound_temp.append((db_cpd_id, 'ATLAS'))

    elif cpd in currentcpds and 'c0' not in currentcpds[cpd]:
        compartment_keys = currentcpds[cpd].keys()
        db_cpd_id = currentcpds[cpd][compartment_keys[0]]
        db_cpd_id = re.sub('_'+compartment_keys[0], '_c0', db_cpd_id)
        reaction_compound_temp.append((rxn, db_cpd_id, is_prod, stoich, 0))
        if db_cpd_id not in dbcpds:
            compound_temp.append((db_cpd_id, 'None', 'c0', cpd))
        model_compound_temp.append((db_cpd_id, 'ATLAS'))

    else:
        if inchidb:
            db_cpd_id = get_inchi_4_cpd(cpd)
            if db_cpd_id:
                reaction_compound_temp.append((rxn, db_cpd_id+'_c0', is_prod, stoich, 0))
                if db_cpd_id+'_c0' not in dbcpds:
                    compound_temp.append((db_cpd_id+'_c0', 'None', 'c0', cpd))
                if (cpd+'_c0', db_cpd_id+'_c0') not in original_db_cpd_current:
                    original_db_cpd_temp.append((cpd+'_c0', db_cpd_id+'_c0'))
                model_compound_temp.append((db_cpd_id+'_c0', 'ATLAS'))
            else:
                reaction_compound_temp.append((rxn, cpd+'_c0', is_prod, stoich, 0))
                if cpd+'_c0' not in dbcpds:
                    compound_temp.append((cpd+'_c0', 'None', 'c0', cpd))
                model_compound_temp.append((cpd+'_c0', 'ATLAS'))
        else:
            reaction_compound_temp.append((rxn, cpd+'_c0', is_prod, stoich, 0))
            if cpd+'_c0' not in dbcpds:
                compound_temp.append((cpd+'_c0', 'None', 'c0', cpd))
                model_compound_temp.append((cpd+'_c0', 'ATLAS'))
    return (model_compound_temp, compound_temp, reaction_compound_temp, original_db_cpd_temp)

def get_inchi_4_cpd(cpd):
    '''Get inchi for a compound if inchidb is True'''
    darray = extract_KEGG_data(KEGG+'get/'+cpd)
    inchivalue = ''
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
                                inchivalue = compounds[0].inchi
                            else:
                                inchivalue = None
                        except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                            inchivalue = None
                    else:
                        inchivalue = None
                except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                    inchivalue = None
    else:
        inchivalue = None
    return inchivalue

def extract_KEGG_data(url):
    '''Extract Kegg db info'''
    try:
        data = urllib2.urlopen(url).read()
        darray = data.split('\n')
        return darray
    except (httplib.BadStatusLine, urllib2.URLError):
        return None

def fill_database(cnx, reactions, reaction_reversible, model_reaction,
                  reaction_protein, reaction_genes, model_compound,
                  compound, reaction_compound, original_db_cpd_new):
    '''fill database with ATLAS information'''
    Q = cnx.execute("SELECT ID FROM model WHERE ID = 'ATLAS'")
    hits = Q.fetchall()
    if not hits:
        print ('STATUS: Adding ATLAS model to model table')
        cnx.execute("INSERT INTO model VALUES (?, ?)", ('ATLAS', 'ATLAS_database'))
        Q = cnx.execute('''SELECT DISTINCT cluster_num FROM cluster''')
        hits = Q.fetchall()
        uniq_clusters = [i[0] for i in hits]
        cnx.execute("INSERT INTO cluster VALUES (?,?)", (len(uniq_clusters)+1, 'ATLAS'))
        cnx.commit()

    cnx.executemany("INSERT INTO model_compound VALUES (?,?)", model_compound)
    cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", model_reaction)
    cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", reaction_protein)
    cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", reaction_genes)
    cnx.executemany("INSERT INTO reaction_reversibility VALUES (?,?)",
                         reaction_reversible)
    cnx.executemany("INSERT INTO reaction VALUES (?,?,?,?)", reactions)
    cnx.executemany("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                         reaction_compound)
    cnx.executemany("INSERT INTO compound VALUES (?,?,?,?)", compound)
    if original_db_cpd_new:
        cnx.executemany("INSERT INTO original_db_cpdIDs VALUES (?,?)",
                             original_db_cpd_new)
    cnx.commit()