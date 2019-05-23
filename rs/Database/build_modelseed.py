import os
import shutil
import re
import sqlite3
import urllib2
import httplib
# from tqdm import tqdm
from multiprocessing import Process, Queue
from requests.exceptions import ConnectionError
from cobra import io
from copy import deepcopy
import pubchempy
import Database
from Database import mackinac
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

KEGG = 'http://rest.kegg.jp/'
PATH = os.path.dirname(os.path.abspath(__file__))
def verbose_print(verbose, line):
    '''verbose print function'''
    if verbose:
        print(line)

def extract_KEGG_data(url, verbose):
    '''Extract Kegg db info'''
    try:
        verbose_print(verbose, url)
        data = urllib2.urlopen(url).read()
        darray = data.split('\n')
        return(darray)
    except urllib2.HTTPError:
        return(None)

def retrieve_exact_inchi_values(new_cpd_keggid, raw_cpd_keggid, cpd_name, compart_info,
                                inchi_pubchem, inchi_cf, inchi_cas, CT, INCHI, verbose):
    '''Retrieve InChI values'''
    if new_cpd_keggid not in inchi_pubchem:
        if raw_cpd_keggid != 'None':
            count = 0
            while count < 3:
                inchi, cas = kegg2pubcheminchi(raw_cpd_keggid, verbose)
                if inchi:
                    count = 4
                else:
                    count+=1
                    verbose_print(verbose, 'WARNING: No INCHI was found from {} will try again for {} time'.format(raw_cpd_keggid, count))
            if not inchi:
                verbose_print(verbose, 'STATUS: trying to inchi for compound {} through pubchem'.format(raw_cpd_keggid))
                inchi, iupac_name, cas = CT.translate(cpd_name)
        else:
            inchi, iupac_name, cas = CT.translate(cpd_name)

        if inchi is not None:
            mol = INCHI.loadMolecule(inchi)
            cf = mol.grossFormula()
            cf = re.sub(' ', '', cf)
            # fp = mol.fingerprint('full')
            # buffer = fp.toBuffer()
            # buffer_array = [str(i) for i in buffer]
            # buffer_string = ','.join(buffer_array)
            inchi_cf[new_cpd_keggid] = cf
            inchi_cas[new_cpd_keggid] = cas
            inchi_pubchem[new_cpd_keggid] = inchi+'_'+compart_info+'0'
        else:
            inchi_cf[new_cpd_keggid] = 'None'
            inchi_cas[new_cpd_keggid] = 'None'
            # inchi_fp[new_cpd_keggid] = 'None'
            inchi_pubchem[new_cpd_keggid] = new_cpd_keggid
    return (inchi_pubchem, inchi_cf, inchi_cas)

def kegg2pubcheminchi(cpd, verbose):
    '''Convvert kegg ID to InChI value'''
    darray = extract_KEGG_data(KEGG+'get/'+cpd, verbose)
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
                    verbose_print(verbose, 'WARNING: Could not get substance for {} {}'.format(sid, cpd))
                    pass
            if 'CAS:' in array:
                index = array.index('CAS:')
                cas = array[index+1]

            
    return (inchicpd, cas)

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

def get_KEGG_IDs(ID, compartment, KEGGdict):
    '''Retrieve KEGG IDs'''
    new_ID = re.sub('_'+compartment+'$', '', ID)
    try:
        KEGG_ID = KEGGdict[new_ID]
        original_KEGG_ID = str(deepcopy(KEGG_ID))
        if not KEGG_ID:
            KEGG_ID = ID+'0'
        else:
            KEGG_ID = str(KEGG_ID)+'_'+str(compartment)+'0'
    except KeyError:
        KEGG_ID = ID+'0'
        original_KEGG_ID = 'None'
    return(KEGG_ID, original_KEGG_ID)

def build_patric_models(genome_id, genome_name, media, username):
    '''Builds patric models on the patric server'''
    try:
        patric_model=mackinac.create_patric_model(str(genome_id), genome_name,
                                                  media_reference='/chenry/public/modelsupport/media/{}'.format(media)) #Carbon-D-Glucose'
    except ConnectionError:
        pass
    except ValueError:
        print ('STATUS: No model identified for genome {} {}'.format(genome_id, genome_name))
    except:
        print ('STATUS: Unexpected error generating patric model for genome {} {}'.format(genome_id, genome_name))

def generate_sbml_output_folder(output_folder):
    '''Generate output folder for sbml fba models if user has specified this option'''
    try:
        os.mkdir(output_folder)
    except OSError:
        shutil.rmtree(output_folder)
        os.mkdir(output_folder)   

class BuildModelSeed(object):
    def __init__(self, username, password,  rxntype, inchidb, DBpath, output_folder,
                 media='Compelete', newdb=True,
                 tokentype='patric', sbml_output=False,
                 processors=4, verbose=False, 
                 patricfile=PATH+'/data/PATRIC_genome_complete_07152018.csv',
                 previously_built_patric_models=False):
        '''initialize'''
        self.username = username
        self.password = password
        self.inchidb = inchidb
        self.tokentype = tokentype
        self.rxntype = rxntype
        self.sbml_output = sbml_output
        self.output_folder = output_folder
        self.media = media
        self.newdb = newdb
        self.processors = processors
        self.verbose = verbose
        self.patricfile = patricfile
        self.previously_built_patric_models = previously_built_patric_models
        self.CPD2KEGG = open_translation_file(PATH+'/data/KbasetoKEGGCPD.txt')
        self.RXN2KEGG = open_translation_file(PATH+'/data/KbasetoKEGGRXN.txt')
        self.inchi_dict = {}
        self.inchi_fp_dict = {}
        self.inchi_cf_dict = {}
        self.inchi_cas_dict = {}
        self.allcompounds = set()
        self.allcompounds_check = set()
        self.allreactions = set()
        self.allreactions_id = set()
        self.reaction_reversibility = {}
        self.CT = pit.CompoundTranslator()
        self.L2D = LoadIntoDB(DBpath, self.verbose, self.inchidb)
        self.originalIDs = set()
        if self.inchidb:
            self.IN = indigo.Indigo()
            self.INCHI = indigo_inchi.IndigoInchi(self.IN)
        token = mackinac.get_token(username=self.username, password=self.password, token_type=self.tokentype)
        if self.sbml_output:
           generate_sbml_output_folder(self.output_folder+'sbml_models/')
        self.load_complete_genomes()
        self.get_model_from_patric()

    def process_cobra_model(self, model, genome_id):
        '''Adds information from cobra model into appropriate arrays'''
        print ('STATUS: Processing cobra model {}'.format(genome_id))
        model_compounds = []
        model_reactions = []
        reaction_gene = []
        reaction_protein = []
        reactions_compounds = []
        cpd_IDs = {}
        print ('STATUS: Getting metabolites for cobra model {}'.format(genome_id))
        for cpd in model.metabolites:
            new_cpd_keggid, raw_cpd_keggid=get_KEGG_IDs(cpd.id, cpd.compartment, self.CPD2KEGG)
            cpd_IDs[cpd.id] = new_cpd_keggid
            if self.inchidb:
                self.inchi_dict, self.inchi_cf_dict,self.inchi_cas_dict = retrieve_exact_inchi_values(new_cpd_keggid, raw_cpd_keggid, cpd.name,
                                                                                                      cpd.compartment, self.inchi_dict,
                                                                                                      self.inchi_cf_dict, self.inchi_cas_dict, 
                                                                                                      self.CT, self.INCHI, self.verbose)
                model_compounds.append((self.inchi_dict[new_cpd_keggid], genome_id))
                if self.inchi_dict[new_cpd_keggid] not in self.allcompounds_check:
                    self.allcompounds.add((self.inchi_dict[new_cpd_keggid], cpd.name, cpd.compartment+'0', raw_cpd_keggid,
                                        self.inchi_cf_dict[new_cpd_keggid], self.inchi_cas_dict[new_cpd_keggid]))
                    if self.inchi_dict[new_cpd_keggid].startswith('InChI'):
                        self.originalIDs.add((new_cpd_keggid, self.inchi_dict[new_cpd_keggid]))
                    self.allcompounds_check.add(self.inchi_dict[new_cpd_keggid])
            else:
                model_compounds.append((new_cpd_keggid, genome_id))
                self.allcompounds.add((new_cpd_keggid, cpd.name, cpd.compartment+'0', raw_cpd_keggid, 'None', 'None'))

        print ('STATUS: Getting reactions for cobra model {}'.format(genome_id))
        for reaction in model.reactions:
            try:
                new_rxn_keggid, raw_rxn_keggid=get_KEGG_IDs(reaction.id, list(reaction.compartments)[0], self.RXN2KEGG)
            except (KeyError, IndexError):
                print ('STATUS: {} has no compartments {} make compartment c0'.format(reaction.id, reaction.compartments))
                new_rxn_keggid, raw_rxn_keggid=get_KEGG_IDs(reaction.id, 'c', self.RXN2KEGG)
            if new_rxn_keggid.startswith('bio'):
                new_rxn_keggid  = new_rxn_keggid+'_'+str(genome_id)
            model_reactions.append((new_rxn_keggid, genome_id, str(reaction.reversibility).lower()))
            for gene in reaction.genes:
                reaction_gene.append((str(new_rxn_keggid), str(genome_id), str(gene)))
            ECnames = reaction.gene_name_reaction_rule
            match = re.findall('\(EC \d+\.\d+\.\d+\.\d+\)', ECnames)
            match_set = set(match)
            if match_set:
                for ec in match_set:
                    reaction_protein.append((str(new_rxn_keggid), str(genome_id), str(ec)))
            if new_rxn_keggid not in self.allreactions_id:
                for cpd in reaction.reactants:
                    stoich = reaction.get_coefficient(cpd.id)
                    if self.inchidb:
                        reactions_compounds.append((new_rxn_keggid, self.inchi_dict[cpd_IDs[cpd.id]], 0, abs(stoich), 0))
                    else:
                        reactions_compounds.append((new_rxn_keggid, cpd_IDs[cpd.id], 0, abs(stoich), 0))
                for cpd in reaction.products:
                    stoich = reaction.get_coefficient(cpd.id)
                    if self.inchidb:
                        reactions_compounds.append((new_rxn_keggid, self.inchi_dict[cpd_IDs[cpd.id]], 1, abs(stoich), 0))
                    else:
                        reactions_compounds.append((new_rxn_keggid, cpd_IDs[cpd.id], 1, abs(stoich), 0))
                self.allreactions.add((new_rxn_keggid, reaction.name, raw_rxn_keggid, self.rxntype))
                self.allreactions_id.add(new_rxn_keggid)
            if new_rxn_keggid not in self.reaction_reversibility:
                self.reaction_reversibility[new_rxn_keggid] = str(reaction.reversibility).lower()
            else:
                if self.reaction_reversibility[new_rxn_keggid] == 'false' and str(reaction.reversibility).lower() == 'true':
                    self.reaction_reversibility[new_rxn_keggid] = str(reaction.reversibility).lower()
        #ADD TO DATABASE

        self.L2D.add_model_compounds(model_compounds)
        self.L2D.add_model_reactions(model_reactions, reaction_gene, reaction_protein)
        self.L2D.add_reaction_compound(reactions_compounds, self.newdb)

    def load_complete_genomes(self):
        '''Loads a list of patric genome IDs that are complete genomes'''
        self.complete_genomes = []
        print ('STATUS: Loading complete patric genome IDs')
        with open(self.patricfile) as fin:
            header = fin.readline()
            for line in fin:
               larray = line.split(',') 
               self.complete_genomes.append((larray[0], larray[1]))

    def get_model_from_patric(self):
        '''Builds models in patric and converts them to cobra models which are then imported into retsynth database'''
        model_compartments = set()
        model_ids = []
        
        print ('STATUS: getting patric modeling information from genomes')
        models_in_ws = mackinac.list_patric_models()
        if models_in_ws is not None:
            models_in_ws_id = [m['id'] for m in models_in_ws]
        else:
            models_in_ws_id = []
        new_patric_models = []
        count_models = 0
        self.complete_genomes_new_name = []
        for (genome_id, genome_name) in self.complete_genomes:
            genome_name = re.sub("""\'|\"""", '', genome_name)
            genome_name = re.sub(' ', '_', genome_name)
            self.complete_genomes_new_name.append((genome_id, genome_name+'_'+self.media))
            print ("STATUS: getting patric model for genome {} {}".format(genome_name+'_'+self.media, genome_id))
            if genome_name+'_'+self.media not in models_in_ws_id:
                new_patric_models.append((genome_id, genome_name+'_'+self.media, self.media, self.username))
            else:
                count_models+=1
        print ('STATUS: {} models already generated'.format(count_models))
        ###Multiple processors used to get patric models 
        if not self.previously_built_patric_models:
            print ("STATUS: Generating new patric models")
            model_chunks = [new_patric_models[i:i+self.processors] for i in range(0, len(new_patric_models), self.processors)]
            for models in model_chunks:
                processes = []
                for model in models:
                    processes.append(Process(target=build_patric_models, args=(model[0], model[1], model[2], model[3])))
                for p in processes:
                    p.start()
                for p in processes:
                    p.join()
        else:
            print ('STATUS: Using models already in patric')    
        print ('STATUS: Converting patric models into cobra models and loading into database ... ')
        count = 0
        models = set()
        for (genome_id, genome_name) in self.complete_genomes_new_name:
            if genome_id not in models:
                count_attemps = 0
                while count_attemps < 3:
                    try:
                        print ('STATUS: retrieving cobra model for {}'.format(genome_name))
                        cobra_model = mackinac.create_cobra_model_from_patric_model(genome_name)
                        count_attemps = 4
                        models.add(genome_id)
                        self.process_cobra_model(cobra_model, genome_id)
                        model_ids.append((genome_id, genome_name))
                        compartments = cobra_model.compartments
                        for comp_id, name in compartments.iteritems():
                            model_compartments.add((comp_id+'0', name))
                        if self.sbml_output:
                            io.write_sbml_model(cobra_model, self.output_folder+'sbml_models/'+genome_name+'_'+self.media)
                    except ConnectionError:
                        count_attemps+=1
                        print ('STATUS: Unable to get cobra model for genome {} {}, error 1, attempt {}'.format(genome_id, genome_name, count_attemps))
                        pass
                    except mackinac.SeedClient.ObjectNotFoundError:
                        count_attemps+=1
                        print ('STATUS: Unable to get cobra model for genome {} {}, error 2, attempt {}'.format(genome_id, genome_name, count_attemps))
                        pass
                    except Database.mackinac.SeedClient.ObjectNotFoundError:
                        count_attemps+=1
                        print ('STATUS: Unable to get cobra model for genome {} {},  error 3, attempt {}'.format(genome_id, genome_name, count_attemps))
                        pass
                    except:
                        count_attemps+=1
                        print ('STATUS: Unexpected error for cobra model for genome {} {},  error 4, attempt {}'. format(genome_id, genome_name, count_attemps))
            else:
                print ('STATUS: {} genome already in database'.format(genome_id))
            reaction_revers_total = []
        for rxnid, revers in self.reaction_reversibility.iteritems():
            reaction_revers_total.append((rxnid, revers))
        model_ids = list(set(model_ids))
        if self.newdb:
            self.L2D.add_all_info_new(list(self.allcompounds), list(self.originalIDs),
                                      list(self.allreactions),
                                      reaction_revers_total, model_ids,
                                      list(model_compartments))
            self.L2D.add_cluster_info()
        else:
            self.L2D.add_all_info_existing(list(self.allcompounds), list(self.originalIDs),
                                           list(self.allreactions),
                                           reaction_revers_total, model_ids,
                                           list(model_compartments))



class LoadIntoDB(object):
    def __init__(self, DBpath, verbose, inchidb):
        '''initialize'''
        self.verbose = verbose
        self.inchidb = inchidb
        self.cnx = sqlite3.connect(DBpath)
        self.cnx.execute("PRAGMA synchronous = OFF")
        self.cnx.execute("PRAGMA journal_mode = OFF")
        self.cnx.commit()
    
    def add_model_compounds(self, model_compounds):
        '''Loads model compound information into database'''
        self.cnx.executemany("INSERT INTO model_compound VALUES (?,?)", model_compounds)
        self.cnx.commit()

    def add_model_reactions(self, model_reactions, reaction_genes, reaction_protein):
        '''Loads reaction information into database'''
        self.cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", model_reactions)
        self.cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", reaction_genes)
        self.cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", reaction_protein)
        self.cnx.commit()
    
    def add_reaction_compound(self, reaction_compound, newdb):
        '''Loads reaction compound information into database'''
        if newdb:
            self.cnx.executemany("INSERT INTO reaction_compound VALUES (?,?,?,?,?)", reaction_compound)
            self.cnx.commit()
        else:
            q=self.cnx.execute("SELECT ID FROM reaction") 
            hits = q.fetchall()
            db_rxns = [i[0] for i in hits]
            for rxn in reaction_compound:
                if rxn[0] not in db_rxns:
                    self.cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)", rxn)
            self.cnx.commit()          

    def add_all_info_new(self, allcompounds, originalIDs, allreactions, reaction_reversibility, model_ids, model_compartments):
        '''Loads unique compound, reaction, compartment and model information into new database'''
        self.cnx.executemany("INSERT INTO compound VALUES (?,?,?,?,?,?)", allcompounds)
        self.cnx.executemany("INSERT INTO reaction VALUES (?,?,?,?)", allreactions)
        self.cnx.executemany("INSERT INTO reaction_reversibility VALUES (?,?)", reaction_reversibility)
        self.cnx.executemany("INSERT INTO compartments VALUES (?,?)", model_compartments)
        self.cnx.executemany("INSERT INTO model VALUES (?,?)", model_ids)
        self.cnx.executemany("INSERT INTO fba_models VALUES (?,?)", model_ids)
        if self.inchidb:
            self.cnx.executemany("INSERT INTO original_db_cpdIDs VALUES (?,?)", originalIDs)
        self.cnx.commit()

    def add_all_info_existing(self, allcompounds, originalIDs, allreactions, reaction_reversibility, model_ids, model_compartments):
        '''Loads unique compound, reaction, compartment and model information into preexisting database'''
        q=self.cnx.execute("SELECT ID FROM compound") 
        hits = q.fetchall()
        db_cpds = [i[0] for i in hits]
        for cpd in allcompounds:
            if cpd[0] not in db_cpds:
                self.cnx.execute("INSERT INTO compound VALUES (?,?,?,?,?,?)", cpd)
        self.cnx.commit()
        
        q=self.cnx.execute("SELECT ID FROM original_db_cpdIDs") 
        hits = q.fetchall()
        db_cpds = [i[0] for i in hits]
        for cpd in originalIDs:
            if cpd[0] not in db_cpds:
                self.cnx.execute("INSERT INTO original_db_cpdIDs VALUES (?,?)", cpd)
        self.cnx.commit()       
        
        q=self.cnx.execute("SELECT ID FROM reaction") 
        hits = q.fetchall()
        db_rxns = [i[0] for i in hits]
        for rxn in allreactions:
            if rxn[0] not in db_rxns:
                self.cnx.execute("INSERT INTO reaction VALUES (?,?,?,?)", rxn)
        self.cnx.commit()

        q=self.cnx.execute("SELECT * FROM reaction_reversibility") 
        hits = q.fetchall()
        db_rxns_rv = [i[0] for i in hits]
        db_rxns_rv1 = [i[1] for i in hits]
        for rxn in reaction_reversibility:
            if rxn[0] not in db_rxns_rv:
                self.cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)", rxn)
            else:
                if rxn[1] != db_rxns_rv1[db_rxns_rv.index(rxn[0])] and rxn[0] == 'true':
                    self.cnx.execute("""UPDATE reaction_reversibility SET is_reversible = 'true' WHERE reaction_ID = ?""",(rxn[0],))
        self.cnx.commit()

        q=self.cnx.execute("SELECT ID FROM compartments") 
        hits = q.fetchall()
        db_compartments = [i[0] for i in hits]
        for compartment in model_compartments:
            if compartment[0] not in db_compartments:
                self.cnx.execute("INSERT INTO compartments VALUES (?,?)", compartment)
        self.cnx.commit()

        q=self.cnx.execute("SELECT ID FROM model") 
        hits = q.fetchall()
        db_models = [i[0] for i in hits]
        for model in model_ids:
            if model[0] not in db_models:
                self.cnx.execute("INSERT INTO model VALUES (?,?)", model)
        self.cnx.commit()
        
    def add_cluster_info(self):
        '''Loads cluster information into database'''
        print ('STATUS: Adding patric cluster info')
        q=self.cnx.execute("SELECT * FROM cluster")
        db_clusters = q.fetchall()
        metabolic_clusters = {}
        cluster_org = {}
        
        if len(db_clusters) > 0:
            verbose_print(self.verbose, 'STATUS: getting cluster information which is already in the database...')
            for cluster in db_clusters:
                cluster_org.setdefault(cluster[0], []).append(cluster[1])
                if cluster[0] not in metabolic_clusters:
                    metabolic_clusters[cluster[0]]={}
                    sorted_cpds = self.get_model_sorted_cpds(cluster[1])
                    sorted_rxns = self.get_model_sorted_rxns(cluster[1])
                    metabolic_clusters[1]['comps'] = sorted_cpds
                    metabolic_clusters[1]['rxns'] = sorted_rxns                    

        q=self.cnx.execute("SELECT ID FROM model")
        hits = q.fetchall()
        db_models = [i[0] for i in hits]
        for model in db_models:
            sorted_cpds = self.get_model_sorted_cpds(model)
            sorted_rxns = self.get_model_sorted_rxns(model)

            if len(metabolic_clusters) == 0:
                metabolic_clusters[1]={}
                metabolic_clusters[1]['comps'] = sorted_cpds
                metabolic_clusters[1]['rxns'] = sorted_rxns
                cluster_org[1] = []
                cluster_org[1].append(model)
            else:
                count_keys = 0
                for key in metabolic_clusters:
                    if (sorted_cpds == metabolic_clusters[key]['comps'] and
                            sorted_rxns == metabolic_clusters[key]['rxns']):
                        cluster_org[key].append(model)
                        break
                    else:
                        count_keys += 1
                if len(metabolic_clusters) == count_keys:
                    num_keys = len(metabolic_clusters)
                    metabolic_clusters[num_keys+1] = {}
                    metabolic_clusters[num_keys+1]['comps'] = sorted_cpds
                    metabolic_clusters[num_keys+1]['rxns'] = sorted_rxns
                    cluster_org[num_keys+1] = []
                    cluster_org[num_keys+1].append(model)
        for key, value in cluster_org.iteritems():
            for v in value: 
                q = self.cnx.execute("SELECT cluster_num FROM cluster WHERE cluster_num= ? AND ID=?",
                                (key, v))
                result = q.fetchone()
                if result is None:     
                    self.cnx.execute("INSERT INTO cluster VALUES (?,?)", (key, v))
            self.cnx.commit()

    def get_model_sorted_cpds(self, model):
        q = self.cnx.execute("SELECT cpd_ID FROM model_compound WHERE model_ID=?",(model,))
        hits = q.fetchall()
        cpds = [i[0] for i in hits]
        sorted_cpds = sorted(cpds)
        return(sorted_cpds)

    def get_model_sorted_rxns(self, model):
        q = self.cnx.execute("SELECT reaction_ID FROM model_reaction WHERE model_ID=?",(model,))
        hits = q.fetchall()
        rxns = [i[0] for i in hits]
        sorted_rxns = sorted(rxns)
        return(sorted_rxns)
