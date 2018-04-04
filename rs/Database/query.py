from __future__ import print_function
__author__ = 'Leanne Whitmore and Lucy Chian'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Interface for database of FBA models and compounds within SBML files'

import sqlite3

def test_db_4_error(cnx, query, db):
    try:
        Q = cnx.execute(query)
        return(Q, cnx)
    except sqlite3.DatabaseError:
        print ('WARNING: Database Error with ('+str(query)+') ...reconnect to the database')
        conn = sqlite3.connect(db, check_same_thread=False)
        conn.text_factory = str
        cnx = conn.cursor()
        try:
            Q = cnx.execute(query)
            return(Q, cnx)
        except sqlite3.DatabaseError:
            print ('WARNING: Database Error with ('+str(query)+') ...reconnect to the database four second time')
            conn = sqlite3.connect(db, check_same_thread=False)
            conn.text_factory = str
            cnx = conn.cursor()
            try:
                Q = cnx.execute(query)
                return(Q, cnx)
            except sqlite3.DatabaseError:
                print ('WARNING: Database Error with ('+str(query)+') ...reconnect to the database four third time returning empty query')
                return (None, cnx)

class Connector(object):
    """Connects to a generated database"""
    def __init__(self, db):
        self.db = db
        self.conn = sqlite3.connect(db, check_same_thread=False)
        self.conn.text_factory = str
        self.cnx = self.conn.cursor()

    def get_uniq_metabolic_clusters(self):
        '''Retrieves unique metabolic clusters (organisms
            with the exact same metabolism) in the database'''
        query = "SELECT DISTINCT cluster_num FROM cluster"
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        #Q = self.cnx.execute(query,)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                uniq_clusters = [i[0] for i in hits]
                return uniq_clusters
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_models_from_cluster(self, cluster):
        '''Retrieves model IDs from a specified cluster in the database'''
        query = "select ID from cluster where cluster_num = '%s'" % cluster
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                if hits is not None:
                    results = [i[0] for i in hits]
                    return results
                else:
                    return str(None)
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_all_models(self):
        '''Retrieves all model IDs in the database'''
        query = "select * from model"
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchall()
                if result is not None:
                    return result
                else:
                    return str(None)
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_organism_name(self, organism_ID):
        '''Retrieves name of metabolic model given a specific model ID'''
        query = "select file_name from model where ID = '%s'" % organism_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    return str(result[0])
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_organism_ID(self, organism_name):
        '''Retrieves ID of metabolic model given a specific model name'''
        organism_name = str('%')+organism_name+str('%')
        query = "select ID from model where name like '%s'" % organism_name
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    return str(result[0])
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_compound_ID(self, compound_ID):
        '''Retrieves compound ID given a compound name'''
        compound_ID = compound_ID+str('%')
        query = "select ID from compound where name like '%s'" % compound_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchall()
                if result is None or len(result) == 0:
                    return str(None)
                else:
                    return result
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_compound_name(self, compound_ID):
        '''Retrieves compound name given a compound ID'''
        if compound_ID.strip() == "":
            return None

        query = "select name from compound where ID = '%s'" % compound_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    try:
                        return str(result[0])
                    except UnicodeEncodeError:
                        return str(result[0].decode('utf-8'))
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_compound_compartment(self, compound_ID):
        '''Retrieves the compartment that the compound is in'''
        if compound_ID.strip() == "":
            return None

        query = "select compartment from compound where ID = '%s'" % compound_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    return str(result[0])
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_reaction_name(self, reaction_ID):
        '''Retrieves name of the reaction given the reaction ID'''
        if reaction_ID.strip() == "":
            return None

        query = "select name from reaction where ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    return str(result[0])
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_reactions(self, compound_ID, is_prod):
        '''Retrieves reaction IDs that have a given compound ID as a reactant or product'''
        if compound_ID.strip() == "":
            return False

        query = "select reaction_ID from reaction_compound where cpd_ID = '%s' and is_prod = '%s'" % (compound_ID, is_prod)
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return False
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_reaction_species(self, reaction_ID):
        ''' Retrieves compound IDs that are in a given a reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select model_ID from model_reaction indexed by\
                 modelreaction_ind2 where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_reactants(self, reaction_ID):
        '''Retrieves reactants (compound IDs) of a given reaction'''
        if reaction_ID.strip() == "":
            return False

        query = "select cpd_ID from reaction_compound where reaction_ID = '%s' and\
                 is_prod = (0)" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return []
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_reactants_reactions(self, compound_ID):
        ''' Retrieves reactions that have a given compound (ID) as a reactant'''
        if compound_ID.strip() == "":
            return False

        query = ("select reaction_ID from reaction_compound where cpd_ID = '%s' and \
                is_prod = (0)") % compound_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return []
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_products_reactions(self, compound_ID):
        ''' Retrieves reactions that have a given compound (ID) as a product'''
        if compound_ID.strip() == "":
            return False

        query = ("select reaction_ID from reaction_compound where cpd_ID = '%s' and \
                is_prod = (1)") % compound_ID

        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return []
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_products(self, reaction_ID):
        '''Retrieves products (compound IDs) of a given reaction'''
        if reaction_ID.strip() == "":
            return False

        query = ("select cpd_ID from reaction_compound where reaction_ID = '%s' and \
                is_prod = (1)") % reaction_ID

        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return []
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_compounds_in_model(self, organism_ID):
        '''Retrives all compounds in a metabolic model given model ID'''
        if organism_ID.strip() == "":
            return False

        query = "select cpd_ID from model_compound indexed by \
                modelcompound_ind1 where model_ID = '%s'" % organism_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return False
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_all_compounds(self):
        '''Retrieves all compounds in the database'''
        query = "select ID from compound"
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return False
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_reactions_in_model(self, organism_ID):
        '''Retrieves all reactions in a metabolic model given model ID'''
        if organism_ID.strip() == "":
            return False

        query = "select reaction_ID from model_reaction indexed by \
                modelreaction_ind1 where model_ID = '%s'" % organism_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return False
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_all_reactions(self):
        '''Retrieves all reactions in the database'''
        query = "select ID from reaction"
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return False
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def is_reversible(self, organism_ID, reaction_ID):
        '''Retrieves reverisbility information of a reaction
            in a specified metabolic model (model ID)'''
        if organism_ID.strip() == "" or reaction_ID.strip() == "":
            return None

        query = "select is_rev from model_reaction where \
                model_ID = '%s' and reaction_ID = '%s'" % (organism_ID, reaction_ID)
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()[0]
                return str(result)
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def is_reversible_all(self, reaction_ID):
        '''Retrieves reversibility information of a reaction independent of model'''
        if reaction_ID.strip() == "":
            return None

        query = "select is_reversible from reaction_reversibility \
                where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchall()[0][0]
                return str(result)
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_genes(self, reaction_ID, organism_ID):
        '''Retrieves gene associations for a reaction
             of a given metabolic network (model ID)'''
        if organism_ID.strip() == "" or reaction_ID.strip() == "":
            return None
        query = "select gene_ID from reaction_gene where \
                reaction_ID = '%s' and model_ID = '%s'" % (reaction_ID, organism_ID)
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    return str(result[0])
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_proteins(self, reaction_ID, organism_ID):
        '''Retrieves protein associations for a reaction of a given metabolic network (model ID)'''
        if organism_ID.strip() == "" or reaction_ID.strip() == "":
            return None
        query = "select protein_ID from reaction_protein where\
                 reaction_ID = '%s' and model_ID = '%s'" % (reaction_ID, organism_ID)
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is None:
                    return str(None)
                else:
                    return str(result[0])
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_stoichiometry(self, reaction_ID, compound_ID, is_prod):
        '''Retrieves stoichiometry of a compound for a given reaction'''
        if compound_ID.strip() == "" or reaction_ID.strip() == "":
            return None
        query = "select stoichiometry from reaction_compound where \
                reaction_ID = '%s' and cpd_ID = '%s' and is_prod = '%s'" % (reaction_ID, compound_ID, is_prod)
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                return result
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_catalysts(self, reaction_ID):
        '''Retrieves the catalyst of reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select catalysts_ID from reaction_catalysts where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_compartment(self, compartment):
        '''Retrieves the compartment ID'''
        if compartment.strip() == "":
            return None
        compartment = str('%')+compartment+str('%')
        query = "select ID from compartments where name like '%s'" % compartment
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return(str(None))

    def get_solvents(self, reaction_ID):
        '''Retrieves solvent of reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select solvents_ID from reaction_solvents where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_temperature(self, reaction_ID):
        '''Retrieves the temperaature reaction is performed at'''
        if reaction_ID.strip() == "":
            return None

        query = "select temperature from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_pressure(self, reaction_ID):
        '''Retrieves the pressure reaction is performed at'''
        if reaction_ID.strip() == "":
            return None

        query = "select pressure from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_time(self, reaction_ID):
        '''Retrieves the time that is required to perform reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select total_time from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_yield(self, reaction_ID):
        '''Retrieves yield that was reported with reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select yield from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_reference(self, reaction_ID):
        '''Retrieves the reference of reaction'''
        if reaction_ID.strip() == "":
            return None

        query = "select reference from reaction_spresi_info where reaction_ID = '%s'" % reaction_ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'
        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_reactions_based_on_type(self, rxntype):
        '''Retrieves reactions based on type'''
        if rxntype.strip() == "":
            return None

        query = "select ID from reaction where type = '%s'" % rxntype
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)

        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'

        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_all_keggIDs(self):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from reaction" 
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return 'None'

        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_kegg_reaction_ID(self, ID):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from reaction where ID = '%s'"  % ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is not None:
                    return result[0]
                else:
                    return 'None'
            else:
                return 'None'

        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_kegg_cpd_ID(self, ID):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from compound where ID = '%s'"  % ID
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                result = Q.fetchone()
                if result is not None:
                    return result[0]
                else:
                    return 'None'
            else:
                return 'None'

        else:
            print ('WARNING: Issue with database')
            return str(None)

    def get_all_kegg_cpd_ID(self):
        '''Retrieves reactions based on type'''
        query = "select kegg_id from compound"
        Q, self.cnx = test_db_4_error(self.cnx, query, self.db)
        if Q:
            if Q.arraysize > 0:
                hits = Q.fetchall()
                return [i[0] for i in hits]
            else:
                return str(None)
        else:
            print ('WARNING: Issue with database')
            return str(None)