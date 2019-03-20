from __future__ import print_function
__author__ = 'Leanne Whitmore and Lucy Chian'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Interface for database of FBA models and compounds within SBML files'

import sqlite3
import time

def test_db_4_error(conn, cnx, query, db, count):
    count+=1
    while count < 15:
        try:
            Q = cnx.execute(query)
            return(Q, cnx)
        except sqlite3.DatabaseError:
            print ('WARNING: Database Error with ('+str(query)+') ...reconnect to the database for {} time'.format(count))
            conn.close()
            conn = sqlite3.connect(db, check_same_thread=False)
            conn.text_factory = str
            cnx = conn.cursor()
            # cnx.execute("PRAGMA synchronous = NORMAL")
            cnx.execute("PRAGMA journal_mode = OFF")
            conn.commit()
            Q, cnx = test_db_4_error(conn, cnx, query, db, count)
            return(Q, cnx)
    if count >= 15:
        print ('WARNING: Database could not get query {} therefore returning None at {} time'.format(query, count))
        return ('Errored', cnx)

def fetching_all_query_results(Q, conn, cnx, db, query, count):
    count+=1
    while count < 15:
        if Q is not None and Q != 'Errored':
            try:
                hits = Q.fetchall()
                return hits
            except sqlite3.DatabaseError:
                print ('WARNING: Database Error with fetching all for query {} results...reconnect to the database for {} time'.format(count, query))
                time.sleep(5)
                Q, cnx = test_db_4_error(conn, cnx, query, db, count)
                hits = fetching_all_query_results(Q, conn, cnx, db, query, count)
                return hits
        else:
            Q, cnx = test_db_4_error(conn, cnx, query, db, count)
            hits = fetching_all_query_results(Q, conn, cnx, db, query, count)
            return hits       
    if count >= 15:
        print ('WARNING: could not fetch results for query {}'.format(query))
        return 'Errored'

def fetching_one_query_results(Q, conn, cnx, db, query, count):
    count+=1
    while count < 15:
        if Q is not None and Q != 'Errored':
            try:
                hits = Q.fetchone()
                return hits
            except sqlite3.DatabaseError:
                print ('WARNING: Database Error with fetching one query results...reconnect to the database for {} time'.format(count))
                time.sleep(5)
                Q, cnx = test_db_4_error(conn, cnx, query, db, count)
                hits = fetching_one_query_results(Q, conn, cnx, db, query, count)
                return hits
        else:
            Q, cnx = test_db_4_error(conn, cnx, query, db, count)
            hits = fetching_one_query_results(Q, conn, cnx, db, query, count)
            return hits
    if count >= 15:
        print ('WARNING: could not fetch results for query {} is the {} time therefore returning None'.format(query, count))
        return 'Errored'

class Connector(object):
    """Connects to a generated database"""
    def __init__(self, database):
        self.database = database

    def connect_to_database(self):
        conn = sqlite3.connect(self.database, check_same_thread=False)
        conn.text_factory = str
        cnx = conn.cursor()
        # cnx.execute("PRAGMA synchronous = NORMAL")
        #cnx.execute("PRAGMA journal_mode = OFF")
        conn.commit()
        return conn, cnx
        
    def get_uniq_metabolic_clusters(self):
        '''Retrieves unique metabolic clusters (organisms
            with the exact same metabolism) in the database'''
        query = "SELECT DISTINCT cluster_num FROM cluster"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits:
            return[i[0] for i in hits]
        else:
            return str(None)


    def get_models_from_cluster(self, cluster):
        '''Retrieves model IDs from a specified cluster in the database'''
        query = "select ID from cluster where cluster_num = '%s'" % cluster
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' or hits is not None:
            return ([i[0] for i in hits])
        else:
            return str(None)

    def get_all_models(self):
        '''Retrieves all model IDs in the database'''
        query = "select * from model"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' or hits is not None:
            return [i for i in hits]
        else:
            return str(None)


    def get_organism_name(self, organism_ID):
        '''Retrieves name of metabolic model given a specific model ID'''
        if organism_ID.strip() == "":
            return None
            
        query = "select file_name from model where ID = '%s'" % organism_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' or hits is not None:
            return str(hits[0])
        else:
            return str(None)
 
    def get_organism_ID(self, organism_name):
        '''Retrieves ID of metabolic model given a specific model name'''
        if organism_name.strip() == "":
            return None
 
        organism_name = str('%')+organism_name+str('%')
        query = "select ID from model where name like '%s'" % organism_name
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' or hits is not None:
            return str(hits[0])
        else:
            return str(None)

    def get_compound_ID(self, compound_name, strict=False):
        '''Retrieves compound ID given a compound name'''        
        if strict:
            query = "SELECT ID FROM compound WHERE name = '%s'" % compound_name
        else:
            compound_name = compound_name+str('%')
            query = "select ID from compound where name like '%s'" % compound_name
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None and len(hits)!=0:
            return hits
        else:
            return str(None)

    def get_compound_name(self, compound_ID):
        '''Retrieves compound name given a compound ID'''
        if compound_ID.strip() == "":
            return None

        query = "select name from compound where ID = '%s'" % compound_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:    
            try:
                return str(hits[0])
            except UnicodeEncodeError:
                return str(hits[0].decode('utf-8'))
        else:
            return str(None)

    def get_compound_compartment(self, compound_ID):
        '''Retrieves the compartment that the compound is in'''
        if compound_ID.strip() == "":
            return None

        query = "select compartment from compound where ID = '%s'" % compound_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return str(hits[0])
        else:
            return str(None)

    def get_reaction_name(self, reaction_ID):
        '''Retrieves name of the reaction given the reaction ID'''
        if reaction_ID.strip() == "":
            return None

        query = "select name from reaction where ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
           return str(hits[0]) 
        else:   
            return str(None)

    def get_reactions(self, compound_ID, is_prod):
        '''Retrieves reaction IDs that have a given compound ID as a reactant or product'''
        if compound_ID.strip() == "":
            return False

        query = "select reaction_ID from reaction_compound where cpd_ID = '%s' and is_prod = '%s'" % (compound_ID, is_prod)
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return None


    def get_reaction_species(self, reaction_ID):
        ''' Retrieves compound IDs that are in a given a reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select model_ID from model_reaction indexed by\
                 modelreaction_ind2 where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_reactants(self, reaction_ID):
        '''Retrieves reactants (compound IDs) of a given reaction'''
        if reaction_ID.strip() == "":
            return False

        query = "select cpd_ID from reaction_compound where reaction_ID = '%s' and is_prod = (0)" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored':
            if hits:   
                return [i[0] for i in hits]
            else:
                return []
        else:
            return 'Errored'

    def get_reactants_reactions(self, compound_ID):
        ''' Retrieves reactions that have a given compound (ID) as a reactant'''
        if compound_ID.strip() == "":
            return None

        query = ("select reaction_ID from reaction_compound where cpd_ID = '%s' and is_prod = (0)") % compound_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return None

    def get_products_reactions(self, compound_ID):
        ''' Retrieves reactions that have a given compound (ID) as a product'''
        if compound_ID.strip() == "":
            return None

        query = ("select reaction_ID from reaction_compound where cpd_ID = '%s' and is_prod = (1)") % compound_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:    
            return [i[0] for i in hits]
        else:
            return None

    def get_products(self, reaction_ID):
        '''Retrieves products (compound IDs) of a given reaction'''
        if reaction_ID.strip() == "":
            return None

        query = ("select cpd_ID from reaction_compound where reaction_ID = '%s' and is_prod = (1)") % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored':
            if hits:   
                return [i[0] for i in hits]
            else:
                return []
        else:
            return 'Errored'

    def get_compounds_in_model(self, organism_ID):
        '''Retrives all compounds in a metabolic model given model ID'''
        if organism_ID.strip() == "":
            return None

        query = "select cpd_ID from model_compound indexed by modelcompound_ind1 where model_ID = '%s'" % organism_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return None

    def get_all_compounds(self):
        '''Retrieves all compounds in the database'''
        query = "select ID from compound"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return None

    def get_reactions_in_model(self, organism_ID):
        '''Retrieves all reactions in a metabolic model given model ID'''
        if organism_ID.strip() == "":
            return None

        query = "select reaction_ID from model_reaction indexed by modelreaction_ind1 where model_ID = '%s'" % organism_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return None

    def get_all_reactions(self):
        '''Retrieves all reactions in the database'''
        query = "select ID from reaction"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return None

    def is_reversible(self, organism_ID, reaction_ID):
        '''Retrieves reverisbility information of a reaction
            in a specified metabolic model (model ID)'''
        if organism_ID.strip() == "" or reaction_ID.strip() == "":
            return None

        query = "select is_rev from model_reaction where model_ID = '%s' and reaction_ID = '%s'" % (organism_ID, reaction_ID)
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return str(hits[0])
        else:
            return str(None) 

    def is_reversible_all(self, reaction_ID):
        '''Retrieves reversibility information of a reaction independent of model'''
        if reaction_ID.strip() == "":
            return None

        query = "select is_reversible from reaction_reversibility \
                where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return str(hits[0][0])
        else:
            return str(None)

    def get_genes(self, reaction_ID, organism_ID):
        '''Retrieves gene associations for a reaction
             of a given metabolic network (model ID)'''
        if organism_ID.strip() == "" or reaction_ID.strip() == "":
            return None
        query = "select gene_ID from reaction_gene where reaction_ID = '%s' and model_ID = '%s'" % (reaction_ID, organism_ID)
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return str(hits[0])
        else:
            return str(None)

    def get_proteins(self, reaction_ID, organism_ID):
        '''Retrieves protein associations for a reaction of a given metabolic network (model ID)'''
        if reaction_ID.strip() == "" or organism_ID.strip() == "":
            return None
        # query = "select protein_ID from reaction_protein where\
        #          reaction_ID = '%s' and model_ID = '%s'" % (reaction_ID, organism_ID)
        query = "SELECT protein_ID from reaction_protein WHERE reaction_ID='%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return str(hits[0])         
        else:
          return str(None) 

    def get_stoichiometry(self, reaction_ID, compound_ID, is_prod):
        '''Retrieves stoichiometry of a compound for a given reaction'''
        if compound_ID.strip() == "" or reaction_ID.strip() == "":
            return None
        query = "select stoichiometry from reaction_compound where reaction_ID = '%s' and cpd_ID = '%s' and is_prod = '%s'" % (reaction_ID, compound_ID, is_prod)
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return hits
        else:
            return str(None)

    def get_catalysts(self, reaction_ID):
        '''Retrieves the catalyst of reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select catalysts_ID from reaction_catalysts where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_compartment(self, compartment):
        '''Retrieves the compartment ID'''
        if compartment.strip() == "":
            return None
        compartment = str('%')+compartment+str('%')
        query = "select ID from compartments where name like '%s'" % compartment
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_solvents(self, reaction_ID):
        '''Retrieves solvent of reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select solvents_ID from reaction_solvents where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_temperature(self, reaction_ID):
        '''Retrieves the temperaature reaction is performed at'''
        if reaction_ID.strip() == "":
            return None

        query = "select temperature from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)


    def get_pressure(self, reaction_ID):
        '''Retrieves the pressure reaction is performed at'''
        if reaction_ID.strip() == "":
            return None

        query = "select pressure from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_time(self, reaction_ID):
        '''Retrieves the time that is required to perform reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select total_time from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_yield(self, reaction_ID):
        '''Retrieves yield that was reported with reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select yield from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_reference(self, reaction_ID):
        '''Retrieves the reference of reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select reference from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_reactions_based_on_type(self, rxntype):
        '''Retrieves reactions based on type'''
        if rxntype.strip() == "":
            return None

        query = "select ID from reaction where type = '%s'" % rxntype
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_reaction_type(self, reaction_ID):
        '''Retrieves reaction type based on ID'''
        if reaction_ID.strip() == "":
            return None

        query = "select type from reaction where ID = '%s'" % reaction_ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_all_keggIDs(self):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from reaction" 
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_kegg_reaction_ID(self, ID):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from reaction where ID = '%s'"  % ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return hits[0]
        else:
            return str(None)

    def get_kegg_cpd_ID(self, ID):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from compound where ID = '%s'"  % ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return hits[0]
        else:
            return str(None)

    def get_all_kegg_cpd_ID(self):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from compound"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    # def get_all_cpd_fp(self):
    #     '''Retrieves all fingerprints'''
    #     query = "select fingerprint from compound"
    #     conn, cnx = self.connect_to_database()
    #     Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
    #     hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
    #     if hits != 'Errored' and hits is not None:
    #         return [i[0] for i in hits]
    #     else:
    #         return str(None)
 
    # def get_cpd_fp(self, ID):
    #     '''Retrieves fingerprint for compound ID'''
    #     query = "select fingerprint from compound where ID = '%s'" % ID
    #     conn, cnx = self.connect_to_database()
    #     Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
    #     hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
    #     if hits != 'Errored' and hits is not None:
    #         return hits[0]
    #     else:
    #         return str(None)

    def get_all_cpd_chemicalformulas(self):
        '''Retrieves all chemicalformulas'''
        query = "select chemicalformula from compound"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)
    def get_cpd_chemicalformula(self, ID):
        '''Retrieves chemicalformula for compound ID'''
        query = "select chemicalformula from compound where ID = '%s'" % ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return hits[0]
        else:
            return str(None)

    def get_cpd_casnumber(self, ID):
        '''Retrieves casnumber for compound ID'''
        query = "select casnumber from compound where ID = '%s'" % ID
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return hits[0]
        else:
            return str(None)

    def get_all_cpd_with_chemicalformula(self, cf):
        '''Retrieves chemicalformula for compound ID'''
        query = "select ID from compound where chemicalformula = '%s'" % cf
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return [i[0] for i in hits]
        else:
            return str(None)

    def get_all_cpd_with_search(self, search):
        '''Retrieves compound name for given search term (name/formula)'''
        # query = ("SELECT name FROM compound WHERE name LIKE '%" + search +
        #     "%' OR chemicalformula LIKE '%" + search + "%'")

        query = ("SELECT name, ID FROM compound WHERE " +
                        "name LIKE '%" + search + "%' OR " +
                        "chemicalformula LIKE '%" + search + "%' OR " +
                        "ID LIKE '%" + search + "%' OR " +
                        "casnumber LIKE '%" + search + "%' "
                    "ORDER BY " +
                        # Exact match
                        "CASE WHEN name = '" + search + "' THEN 2 ELSE 0 END + " +
                        "CASE WHEN ID = '" + search + "' THEN 2 ELSE 0 END + " +
                        "CASE WHEN chemicalformula = '" + search + "' THEN 2 ELSE 0 END + " +
                        "CASE WHEN casnumber = '" + search + "' THEN 2 ELSE 0 END + " +
                        # Starts with ...
                        "CASE WHEN name LIKE '" + search + "%' THEN 1 ELSE 0 END + " +
                        "CASE WHEN chemicalformula LIKE '" + search + "%' THEN 1 ELSE 0 END + " +
                        "CASE WHEN ID LIKE '" + search + "%' THEN 1 ELSE 0 END + " +
                        "CASE WHEN casnumber LIKE '" + search + "%' THEN 1 ELSE 0 END " +
                    "DESC"
                )

        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            # return [i[0] for i in hits]
            return hits
        else:
            return str(None)

    def get_model_ID(self, file_name):
        '''Retrieves model ID for given file_name'''
        query = "SELECT ID FROM model WHERE file_name='%s'" % file_name
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_one_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' and hits is not None:
            return hits[0]
        else:
            return str(None)

    def get_all_fba_models(self):
        '''Retrieves all model IDs in the database'''
        query = "select * from fba_models"
        conn, cnx = self.connect_to_database()
        Q, cnx = test_db_4_error(conn, cnx, query, self.database, 0)
        hits = fetching_all_query_results(Q, conn, cnx, self.database, query, 0)
        if hits != 'Errored' or hits is not None:
            return [i for i in hits]
        else:
            return str(None)