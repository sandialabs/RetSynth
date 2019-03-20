from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Build metabolic db using KEGG db'

from multiprocessing import Process, Queue
import re
import sqlite3
import urllib2
import httplib
import pubchempy
from tqdm import tqdm
from Database import query as Q
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

IN = indigo.Indigo()
INCHI = indigo_inchi.IndigoInchi(IN)
KEGG = 'http://rest.kegg.jp/'

def BuildKEGG(types_orgs, inchidb, processors, currentcpds,
             num_organisms='all', num_pathways='all'):
    """Build metabolic database from KEGG DB"""
    types_orgs = types_orgs.lower()
    orgIDs = extract_KEGG_orgIDs(types_orgs, num_organisms)
    
    metabolic_clusters = {}
    pathwayIDs = {}
    pathwayIDs_set = set()
    removeorg = []
    output_queue = Queue()
    args_organism = [orgIDs.keys()[i:i+processors]
                     for i in range(0, len(orgIDs.keys()), processors)]

    print ('STATUS: get KEGG pathway IDs')
    for orgs in tqdm(args_organism):
        processes = []
        for org in orgs:
            processes.append(Process(target=extract_pathwayIDs,
                                     args=(org, num_pathways, output_queue)))
        for p in processes:
            p.start()
        for p in processes:
            result_tuples = output_queue.get()
            pathwayIDs.update(result_tuples[0])
            removeorg.extend(result_tuples[1])

    for orgID in removeorg:
        del orgIDs[orgID]

    cluster_count = 0
    for orgID, pathways in pathwayIDs.iteritems():
        for pathway in pathways:
            pathwayIDs_set.add(pathway)

        if metabolic_clusters:
            check = False
            for key in metabolic_clusters:
                if pathways == pathwayIDs[metabolic_clusters[key][0]]:
                    check = True
                    metabolic_clusters[key].append(orgID)
                    break
            if check is False:
                cluster_count += 1
                metabolic_clusters.setdefault(cluster_count, []).append(orgID)
        else:
            cluster_count += 1
            metabolic_clusters.setdefault(cluster_count, []).append(orgID)

    print ('STATUS: get KEGG reaction IDs')
    output_queue = Queue()
    pathwayIDs_set_list = list(pathwayIDs_set)
    args_pathways = [pathwayIDs_set_list[i:i+processors]
                     for i in range(0, len(pathwayIDs_set_list), processors)]
    pathway_reactionIDs = {}
    noreactions = []
    for pathways in tqdm(args_pathways):
        processes = []
        for pathway in pathways:
            processes.append(Process(target=extract_reactionIDs, args=(pathway, output_queue)))
        for p in processes:
            p.start()
        for p in processes:
            result_tuples = output_queue.get()
            pathway_reactionIDs.update(result_tuples[0])
            noreactions.extend(result_tuples[1])

    reactionIDs = {}
    reactions_set = set()
    for orgID, pathways in pathwayIDs.iteritems():
        reactionIDs[orgID] = []
        for pathway in pathways:
            try:
                for reaction in pathway_reactionIDs[pathway]:
                    reactionIDs[orgID].append(reaction)
                    reactions_set.add(reaction)
            except KeyError:
                print ('WARNING: No reactions for '+str(pathway))

    reactions_set_list = list(reactions_set)
    reactioninfo_final = {}
    inchi_cf = {}
    inchi_cas = {}
    cpd2inchi_final = {}
    compoundinfo_final = {}
    output_queue = Queue()
    args_reactions = [reactions_set_list[i:i+processors]
                      for i in range(0, len(reactions_set_list), processors)]

    print ('STATUS: get KEGG reaction detailed information')
    for reactions in tqdm(args_reactions):
        processes = []
        for reaction in reactions:
            processes.append(Process(target=process_reaction,
                                     args=(reaction, inchidb, compoundinfo_final,
                                           cpd2inchi_final, inchi_cf, inchi_cas,
                                           currentcpds, output_queue)))
        for p in processes:
            p.start()
        for p in processes:
            result_tuples = output_queue.get()
            reactioninfo_final.update(result_tuples[0])
            cpd2inchi_final.update(result_tuples[1])
            inchi_cf.update(result_tuples[2])
            inchi_cas.update(result_tuples[3])
            compoundinfo_final.update(result_tuples[4])

    return(orgIDs, metabolic_clusters, reactioninfo_final,
           reactionIDs, cpd2inchi_final, inchi_cf, inchi_cas,
           compoundinfo_final)

def extract_KEGG_orgIDs(types_orgs, num_organisms):
    """Retrieve organism IDs in KEGG"""
    print ('STATUS: get KEGG '+types_orgs+' organism IDs')
    orgIDs = {}
    darray = extract_KEGG_data(KEGG+'list/organism')
    if types_orgs == 'all':
        for value in tqdm(darray):
            array = value.split('\t')
            try:
               orgIDs[array[1]] = array[2]
            except IndexError:
                pass 
            if num_organisms != 'all':
                if len(orgIDs) == int(num_organisms):
                     break

    else:
        types_orgs = re.sub('\s+', '', types_orgs)
        types_orgs_array = types_orgs.split(',')
        print (types_orgs_array)
        for value in tqdm(darray):
            array = value.split('\t')
            try:
                array[3] = array[3].lower()
                matches = []
                for org in types_orgs_array:
                    match = re.search(org, array[3])
                    if match:
                        matches.append(match)
                if matches:
                    print (array[3])
                    orgIDs[array[1]] = array[2]
            except IndexError:
                pass
            if num_organisms != 'all':
                if len(orgIDs) == int(num_organisms):
                     break
    print ('STATUS: {} organisms from KEGG'.format(len(orgIDs)))
    return orgIDs

def extract_pathwayIDs(orgID, num_pathways, output_queue):
    """Retrieve pathway IDs"""
    pathwayIDs = {}
    removeorg = []
    darray = extract_KEGG_data(KEGG+'list/pathway/'+orgID)
    if darray:
        for count1, value in enumerate(darray):
            array = value.split('\t')
            array[0] = re.sub('path:'+orgID, '', array[0])
            pathwayIDs.setdefault(orgID, []).append(array[0])
            if num_pathways != 'all':
                if count1 == int(num_pathways):
                    break
    else:
        print ('WARNING: Could not get pathways for '+orgID+' so removing it from organisms')
        removeorg.append(orgID)
    output_queue.put((pathwayIDs, removeorg))

def extract_reactionIDs(pathway, output_queue):
    """Extract reactions in pathways"""
    temp_reactionIDs = {}
    temp_noreactions = []
    darray = extract_KEGG_data(KEGG+'link/rn/map'+pathway)
    if darray:
        for value in darray:
            array = value.split('\t')
            try:
                array[1] = re.sub('rn:', '', array[1])
                temp_reactionIDs.setdefault(pathway, []).append(array[1])
            except IndexError:
                pass
    else:
        print ('WARNING: No reactions linked to pathway '+str(pathway))
        temp_noreactions.append(pathway)
    output_queue.put((temp_reactionIDs, temp_noreactions))

def process_reaction(reactionID, inchidb, compoundinfo, cpd2inchi, inchi_cf, inchi_cas, currentcpds, output_queue):
    """Extract reaction info"""
    reactioninfo = {}
    darray = extract_KEGG_data(KEGG+'get/'+reactionID)
    if darray:
        reactioninfo[reactionID] = {}
        reactioninfo[reactionID]['reactants'] = {}
        reactioninfo[reactionID]['products'] = {}
        for value in darray:
            array = value.split()
            if 'NAME' in array:
                reactioninfo[reactionID]['name'] = array[1]
            if 'ENZYME' in array:
                reactioninfo[reactionID]['enzyme'] = array[1]
            if 'EQUATION' in array:
                array.remove('EQUATION')
                equation_string = ' '.join(array)
                if '<=>' in array:
                    reactioninfo[reactionID]['reversible'] = 'true'
                    reaction = equation_string.split(' <=> ')
                    reactants = reaction[0].split(' + ')
                    products = reaction[1].split(' + ')
                    for reactant in reactants:
                        (reactioninfo,
                         cpd2inchi_individual,
                         inchi_cf_individual,
                         inchi_cas_individual,
                         compoundinfo_individual) = process_compound(reactant, reactionID,
                                                                     reactioninfo, False,
                                                                     inchidb, compoundinfo,
                                                                     cpd2inchi, inchi_cf,
                                                                     inchi_cas, currentcpds)
                        cpd2inchi.update(cpd2inchi_individual)
                        inchi_cf.update(inchi_cf_individual)
                        inchi_cas.update(inchi_cas_individual)
                        compoundinfo.update(compoundinfo_individual)
                    for product in products:
                        (reactioninfo,
                         cpd2inchi_individual,
                         inchi_cf_individual,
                         inchi_cas_individual,
                         compoundinfo_individual) = process_compound(product, reactionID,
                                                                     reactioninfo, True,
                                                                     inchidb, compoundinfo,
                                                                     cpd2inchi, inchi_cf,
                                                                     inchi_cas, currentcpds)
                        cpd2inchi.update(cpd2inchi_individual)
                        inchi_cf.update(inchi_cf_individual)
                        inchi_cas.update(inchi_cas_individual)
                        compoundinfo.update(compoundinfo_individual)
                else:
                    reactioninfo[reactionID]['reversible'] = 'false'
                    reaction = equation_string.split(' => ')
                    reactants = reaction[0].split(' + ')
                    products = reaction[1].split(' + ')
                    for reactant in reactants:
                        (reactioninfo,
                         cpd2inchi_individual,
                         inchi_cf_individual,
                         inchi_cas_individual,
                         compoundinfo_individual) = process_compound(reactant, reactionID,
                                                                     reactioninfo, False,
                                                                     inchidb, compoundinfo,
                                                                     cpd2inchi, inchi_cf,
                                                                     inchi_cas, currentcpds)
                        cpd2inchi.update(cpd2inchi_individual)
                        inchi_cf.update(inchi_cf_individual)
                        inchi_cas.update(inchi_cas_individual)
                        compoundinfo.update(compoundinfo_individual)
                    for product in products:
                        (reactioninfo,
                         cpd2inchi_individual,
                         inchi_cf_individual,
                         inchi_cas_individual,
                         compoundinfo_individual) = process_compound(product, reactionID,
                                                                     reactioninfo, True,
                                                                     inchidb, compoundinfo,
                                                                     cpd2inchi, inchi_cf,
                                                                     inchi_cas, currentcpds)
                        cpd2inchi.update(cpd2inchi_individual)
                        inchi_cf.update(inchi_cf_individual)
                        inchi_cas.update(inchi_cas_individual)
                        compoundinfo.update(compoundinfo_individual)
        if 'enzyme' not in reactioninfo[reactionID].keys():
            reactioninfo[reactionID]['enzyme'] = 'None'
        if 'name' not in reactioninfo[reactionID].keys():
            reactioninfo[reactionID]['name'] = 'None'

    output_queue.put((reactioninfo, cpd2inchi, inchi_cf, inchi_cas, compoundinfo))

def process_compound(cpd, reactionID, reactioninfo, is_prod,
                     inchidb, compoundinfo, cpd2inchi, inchi_cf,
                     inchi_cas, currentcpds):
    """Extract compound info"""
    cpd = re.sub('\s*\(\S+\)$', '', cpd)
    cpd = re.sub('^\(\S+\)\s*', '', cpd)

    match = re.search('^\d', cpd)
    match2 = re.search('^\w+\s+', cpd)
    if match:
        stoichiometry = int(match.group(0))
        cpd = re.sub('^\d+ ', '', cpd)
    elif match2:
        stoichiometry = 1
        cpd = re.sub('^\w+\s+', '', cpd)
    else:
        stoichiometry = 1
    if cpd in compoundinfo and not inchidb:
        reactioninfo = add_metabolite(reactionID, cpd, stoichiometry, is_prod, reactioninfo)
    elif cpd not in compoundinfo and not inchidb:
        reactioninfo = add_metabolite(reactionID, cpd, stoichiometry, is_prod, reactioninfo)
        darray = extract_KEGG_data(KEGG+'get/'+cpd)
        count_name = 0
        if darray:
            for value in darray:
                array = value.split()
                if 'NAME' in array:
                    count_name += 1
                    compoundinfo[cpd] = array[1]
            if count_name == 0:
                compoundinfo[cpd] = 'None'
        else:
            compoundinfo[cpd] = 'None'

    elif inchidb and cpd in cpd2inchi.values():
        reactioninfo = add_metabolite(reactionID, cpd2inchi.keys()[cpd2inchi.values().index(cpd)],
                                      stoichiometry, is_prod, reactioninfo)
    else:
        if cpd in currentcpds and 'c0' in currentcpds[cpd]:
            reactioninfo = add_metabolite(reactionID, currentcpds[cpd]['c0'],
                                          stoichiometry, is_prod, reactioninfo)
        else:
            darray = extract_KEGG_data(KEGG+'get/'+cpd)
            if inchidb:
                inchicpd = None
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
                                            mol = INCHI.loadMolecule(inchicpd)
                                            # fp = mol.fingerprint('full')
                                            # buffer = fp.toBuffer()
                                            # buffer_array = [str(i) for i in buffer]
                                            # buffer_string = ','.join(buffer_array)
                                            cf = mol.grossFormula()
                                            cf = re.sub(' ', '', cf)
                                            cpd2inchi[inchicpd] = cpd
                                            # inchi_fp[inchicpd] = buffer_string
                                            inchi_cf[inchicpd] = cf
                                        else:
                                            inchicpd = cpd
                                    except (pubchempy.PubChemHTTPError, httplib.BadStatusLine,
                                            urllib2.URLError):
                                        inchicpd = cpd
                                else:
                                    inchicpd = cpd
                            except (pubchempy.PubChemHTTPError, httplib.BadStatusLine,
                                    urllib2.URLError):
                                print ('WARNING: Could not get substance for '+ sid)
                                inchicpd = cpd
                    if not inchicpd:
                        inchicpd = cpd
                else:
                    inchicpd = cpd
                reactioninfo = add_metabolite(reactionID, inchicpd,
                                              stoichiometry, is_prod, reactioninfo)
            else:
                reactioninfo = add_metabolite(reactionID, cpd,
                                              stoichiometry, is_prod, reactioninfo)
            count_name = 0
            for value in darray:
                array = value.split()
                if 'NAME' in array:
                    count_name += 1
                    if inchidb:
                        compoundinfo[inchicpd] = array[1]
                    else:
                        compoundinfo[cpd] = array[1]
                if 'CAS:' in array:
                    index = array.index('CAS:')
                    cas = array[index+1]
                    if inchidb:
                        inchi_cas[inchicpd] = cas
                    else:
                        inchi_cas[cpd] = cas
            if count_name == 0:
                if inchidb:
                    compoundinfo[inchicpd] = 'None'
                else:
                    compoundinfo[cpd] = 'None'
    return (reactioninfo, cpd2inchi, inchi_cf, inchi_cas, compoundinfo)

def add_metabolite(reactionID, cpd, stoichiometry, is_prod, reactioninfo):
    """add metabolites to dictionary"""
    if is_prod:
        reactioninfo[reactionID]['products'][cpd] = stoichiometry
    else:
        reactioninfo[reactionID]['reactants'][cpd] = stoichiometry
    return reactioninfo

def extract_KEGG_data(url):
    '''Extract Kegg db info'''
    try:
        data = urllib2.urlopen(url).read()
        darray = data.split('\n')
        return darray
    except (httplib.BadStatusLine, urllib2.URLError):
        return None

class CompileKEGGIntoDB(object):
    """Add KEGG info to sqlite database"""
    def __init__(self, database, type_org, inchidb, processors,
                 num_organisms, num_pathways, rxntype, add):
        """Initialize"""
        self.database = database
        self.inchidb = inchidb
        self.rxntype = rxntype
        self.currentcpds = {}
        self.conn = sqlite3.connect(self.database, check_same_thread=False)
        self.conn.text_factory = str
        self.cnx = self.conn.cursor()

        self.DB = Q.Connector(database)
        A = self.cnx.execute('SELECT * from compound')
        results = A.fetchall()
        for i in results:
            self.currentcpds.setdefault(i[3], {})
            self.currentcpds[i[3]].setdefault(i[2], i[0])
        (self.orgIDs,
         self.metabolic_clusters,
         self.reactioninfo,
         self.reactionIDs,
         self.cpd2inchi,
         self.inchi_cf,
         self.inchi_cas,
         self.compoundinfo) = BuildKEGG(type_org, self.inchidb, processors,
                                        self.currentcpds, num_organisms, num_pathways)
        print ('STATUS: Adding kegg info to database')
        if add:
            self.add_to_preexisting_db()
        else:
            self.cnx = sqlite3.connect(database)
            self.fill_new_database()

    def add_to_preexisting_db(self):
        '''Add KEGG info to already developed database'''
        cytosol = '_'+self.DB.get_compartment('cytosol')[0]
        all_current_keggIDs = set(self.DB.get_all_keggIDs())
        all_cytosol_rxn = {}
        QC = self.cnx.execute("""SELECT * from reaction""")
        results = QC.fetchall()
        for result in results:
            if (re.search('c0', result[0]) is not None) and result[2] != 'None':
                all_cytosol_rxn[result[2]] = result[0]

        print ('STATUS: KEGG ..ENTERING MODELS')
        for orgID in self.orgIDs:
            self.cnx.execute("INSERT INTO model VALUES (?,?)", (orgID, self.orgIDs[orgID]+'_KEGG'))
        self.conn.commit()
        print ('STATUS: KEGG ..ENTERING CLUSTERS')
        Q = self.cnx.execute('''SELECT DISTINCT cluster_num FROM cluster''')
        hits = Q.fetchall()
        uniq_clusters = [i[0] for i in hits]
        for key, orgIDs in self.metabolic_clusters.iteritems():
            for orgID in orgIDs:
                self.cnx.execute("INSERT INTO cluster VALUES (?,?)", (len(uniq_clusters)+key, orgID))
        self.conn.commit()

        '''GET CURRENT COMPOUNDS IN DATABASE'''
        Q = self.cnx.execute("""SELECT ID from compound""")
        hits = Q.fetchall()
        dbcpds = [i[0] for i in hits]
        dbcpds = list(set(dbcpds))

        reaction_protein = []
        reaction_gene = []
        model_reactions = []
        model_compounds = []

        print ('STATUS: KEGG ..ENTERING MODEL REACTIONS, GENES, PROTEINS')
        for orgID, reactions in tqdm(self.reactionIDs.iteritems()):
            cpds_all = set()
            for reaction in reactions:
                if reaction in all_cytosol_rxn:
                    try:
                        rxn = self.reactioninfo[reaction]
                        model_reactions.append((all_cytosol_rxn[reaction], orgID, rxn['reversible']))
                        try:
                            reaction_protein.append((all_cytosol_rxn[reaction], orgID, rxn['enzyme']))
                        except KeyError:
                            reaction_protein.append((all_cytosol_rxn[reaction], orgID, 'None'))
                        reaction_gene.append((all_cytosol_rxn[reaction], orgID, 'None'))
                        cpds_all.update(rxn['reactants'])
                        cpds_all.update(rxn['products'])
                    except KeyError:
                        print ('WARNING: No reaction info for {}'.format(reaction))
                else:
                    try:
                        rxn = self.reactioninfo[reaction]
                        model_reactions.append((reaction+cytosol, orgID, rxn['reversible']))
                        try:
                            reaction_protein.append((reaction+cytosol, orgID, rxn['enzyme']))
                        except KeyError:
                            reaction_protein.append((reaction+cytosol, orgID, 'None'))
                        reaction_gene.append((reaction+cytosol, orgID, 'None'))
                        cpds_all.update(rxn['reactants'])
                        cpds_all.update(rxn['products'])
                    except KeyError:
                        print ('WARNING: No reaction info for {}'.format(reaction))                    
            for cpd in cpds_all:
                if cpd.endswith('c0'):
                    model_compounds.append((cpd, orgID))
                else:
                    model_compounds.append((cpd+cytosol, orgID))
        self.cnx.executemany("INSERT INTO model_compound VALUES (?,?)", model_compounds)
        self.cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", model_reactions)
        self.cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", reaction_protein)
        self.cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", reaction_gene)
        self.conn.commit()

        print ('STATUS: KEGG ..ENTERING REACTIONS')
        all_cpds_to_add = set()
        reaction_reversibility = []
        reactions = []
        reaction_compounds = []
        for reaction in self.reactioninfo:
            if reaction not in all_current_keggIDs:
                reaction_reversibility.append((reaction+cytosol, self.reactioninfo[reaction]['reversible']))
                try:
                    reactions.append((reaction+cytosol, self.reactioninfo[reaction]['name'], reaction, self.rxntype))
                except KeyError:
                    reactions.append((reaction+cytosol, 'None', reaction, self.rxntype))
                for reactant in self.reactioninfo[reaction]['reactants']:
                    all_cpds_to_add.add(reactant)
                    if reactant.endswith('c0'):
                        reaction_compounds.append((reaction+cytosol, reactant, 0,
                                                   self.reactioninfo[reaction]['reactants'][reactant], 0))
                    else:
                        reaction_compounds.append((reaction+cytosol, reactant+cytosol, 0,
                                                   self.reactioninfo[reaction]['reactants'][reactant], 0))
                for product in self.reactioninfo[reaction]['products']:
                    all_cpds_to_add.add(product)
                    if product.endswith('c0'):
                        reaction_compounds.append((reaction+cytosol, product, 1,
                                                   self.reactioninfo[reaction]['products'][product], 0))
                    else:
                        reaction_compounds.append((reaction+cytosol, product+cytosol, 1,
                                                   self.reactioninfo[reaction]['products'][product], 0))

        self.cnx.executemany("INSERT INTO reaction_reversibility VALUES (?,?)",
                                reaction_reversibility)
        self.cnx.executemany("INSERT INTO reaction VALUES (?,?,?,?)",
                                reactions)
        self.cnx.executemany("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                reaction_compounds)
        self.conn.commit()
        compounds = []
        print ('STATUS: KEGG ..ENTERING COMPOUNDS')
        for cpd in all_cpds_to_add:
            if cpd.endswith('_c0') or cpd+'_c0' in dbcpds:
                pass
            else:
                if cpd.startswith('InChI'):
                    try:
                        compounds.append((cpd+cytosol, self.compoundinfo[cpd],
                                          self.DB.get_compartment('cytosol')[0],
                                          self.cpd2inchi[cpd],
                                          self.inchi_cf.get(cpd, 'None'),
                                          self.inchi_cas.get(cpd, 'None')))
                    except KeyError:
                        compounds.append((cpd+cytosol, 'None',
                                          self.DB.get_compartment('cytosol')[0],
                                          self.cpd2inchi[cpd],
                                          self.inchi_cf.get(cpd, 'None'),
                                          self.inchi_cas.get(cpd, 'None')))
                else:
                    try:
                        compounds.append((cpd+cytosol, self.compoundinfo[cpd],
                                          self.DB.get_compartment('cytosol')[0],
                                          'None',
                                          self.inchi_cf.get(cpd, 'None'),
                                          self.inchi_cas.get(cpd, 'None')))
                    except KeyError:
                        compounds.append((cpd+cytosol, 'None',
                                          self.DB.get_compartment('cytosol')[0],
                                          'None',
                                          self.inchi_cf.get(cpd, 'None'),
                                          self.inchi_cas.get(cpd, 'None')))
        self.cnx.executemany("INSERT INTO compound VALUES (?,?,?,?,?,?)", compounds)
        self.conn.commit()

    def fill_new_database(self):
        """Fill database"""
        cytosol = '_c0' #define compartment
        self.cnx.execute("INSERT INTO compartments VALUES (?,?)", ('c0', 'cytosol'))
        self.cnx.execute("INSERT INTO compartments VALUES (?,?)", ('e0', 'extracellular'))
        for orgID in self.orgIDs:
            self.cnx.execute("INSERT INTO model VALUES (?,?)", (orgID, self.orgIDs[orgID]+'_KEGG'))
        self.cnx.commit()

        for key, orgIDs in self.metabolic_clusters.iteritems():
            for orgID in orgIDs:
                self.cnx.execute("INSERT INTO cluster VALUES (?,?)", (key, orgID))
        self.cnx.commit()

        for orgID, reactions in self.reactionIDs.iteritems():
            cpds_all = set()
            for reaction in reactions:
                reaction_cm = reaction+cytosol
                self.cnx.execute("INSERT INTO reaction_protein VALUES (?,?,?)",
                                 (reaction_cm, orgID, self.reactioninfo[reaction]['enzyme']))
                self.cnx.execute("INSERT INTO reaction_gene VALUES (?,?,?)",
                                 (reaction_cm, orgID, 'None'))
                self.cnx.execute("INSERT INTO model_reaction VALUES (?,?,?)",
                                 (reaction_cm, orgID, self.reactioninfo[reaction]['reversible']))
                cpds_all.update(self.reactioninfo[reaction]['reactants'])
                cpds_all.update(self.reactioninfo[reaction]['products'])
            for cpd in cpds_all:
                self.cnx.execute("INSERT INTO model_compound VALUES (?,?)", (cpd+cytosol, orgID))
        self.cnx.commit()

        all_cpds_to_add = set()
        for reaction in self.reactioninfo:
            reaction_cm = reaction+cytosol
            self.cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)",
                             (reaction_cm, self.reactioninfo[reaction]['reversible']))
            self.cnx.execute("INSERT INTO reaction VALUES (?,?,?,?)",
                             (reaction_cm, self.reactioninfo[reaction]['name'], reaction, self.rxntype))
            for reactant in self.reactioninfo[reaction]['reactants']:
                self.cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                 (reaction_cm, reactant+cytosol, 0,
                                  self.reactioninfo[reaction]['reactants'][reactant], 0))
                all_cpds_to_add.add(reactant)
            for product in self.reactioninfo[reaction]['products']:
                self.cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                 (reaction_cm, product+cytosol, 1,
                                  self.reactioninfo[reaction]['products'][product], 0))
                all_cpds_to_add.add(product)
        self.cnx.commit()

        if self.inchidb:
            for inchi in self.cpd2inchi:
                self.cnx.execute("INSERT INTO original_db_cpdIDs VALUES (?,?)",
                                 (self.cpd2inchi[inchi]+cytosol, inchi+cytosol))
            self.cnx.commit()

        for cpd in all_cpds_to_add:
            self.cnx.execute("""INSERT INTO compound VALUES (?,?,?,?,?,?)""",
                             (cpd+cytosol, self.compoundinfo[cpd], 'c0', cpd),
                             self.inchi_cf.get(cpd, 'None'), self.inchi_cas.get(cpd, 'None'))
        self.DB = Q.Connector(self.database)
