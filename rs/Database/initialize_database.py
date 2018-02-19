from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Generates Sqlite database'

import sqlite3

class Createdb(object):
    """Generates new sqlite database and creates table and indicies"""
    def __init__(self,database, inchidb):
        self.inchidb = inchidb
        self.database = sqlite3.connect(database)
        self.create_model_tbls()

    def create_model_tbls(self):
        '''
        Creates tables and indicies in SQLite Database
        '''
        try:
            self.database.execute('''CREATE table model
                                            (ID text,file_name text)''')
            self.database.execute('''CREATE table compound
                                            (ID text, name text,compartment text, kegg_id)''')
            self.database.execute('''CREATE table compartments (ID text, name text)''')
            self.database.execute('''CREATE table model_compound
                                            (cpd_ID text, model_ID text)''')
            self.database.execute('''CREATE table reaction
                                            (ID text, name text, kegg_id, type text)''')
            self.database.execute('''CREATE table model_reaction
                                            (reaction_ID text, model_ID text,is_rev bit(1))''')
            self.database.execute('''CREATE table reaction_compound
                                            (reaction_ID text, cpd_ID text, is_prod bit(1), 
                                            stoichiometry int, filenum int)''')
            self.database.execute('''CREATE table reaction_reversibility
                                            (reaction_ID text, is_reversible bit(1))''')
            self.database.execute('''CREATE table reaction_gene
                                            (reaction_ID text, model_ID text, gene_ID text)''')
            self.database.execute('''CREATE table reaction_protein
                                            (reaction_ID text, model_ID text, protein_ID text)''')
            self.database.execute('''CREATE table cluster
                                            (cluster_num text, ID text)''')
            if self.inchidb is True:
                self.database.execute('''CREATE table original_db_cpdIDs
                                                (ID text, inchi_id text)''')
            print ('STATUS: Generating new database ...')
        except sqlite3.OperationalError:
            print ('WARNING: Database already exists, adding new xml files to existing database ...')
        try:
            self.database.execute('''CREATE INDEX reactioncompound_ind1 ON
                                            reaction_compound(reaction_ID,cpd_ID,is_prod)''')
            self.database.execute('''CREATE INDEX reactioncompound_ind2 ON
                                            reaction_compound(cpd_ID,is_prod)''')
            self.database.execute('''CREATE INDEX modelreaction_ind1 ON
                                            model_reaction(model_ID)''')
            self.database.execute('''CREATE INDEX modelreaction_ind2 ON
                                            model_reaction(reaction_ID)''')
            self.database.execute('''CREATE INDEX modelcompound_ind1 ON
                                            model_compound(model_ID)''')
            self.database.execute('''CREATE INDEX modelcompound_ind2 ON
                                            model_compound(cpd_ID)''')
            self.database.execute('''CREATE INDEX model_ind ON
                                            model(ID)''')
            self.database.execute('''CREATE INDEX reaction_ind ON
                                            reaction(ID)''')
            self.database.execute('''CREATE INDEX compound_ind ON
                                            compound(ID)''')
            self.database.execute('''CREATE INDEX reactiongene_ind ON
                                            reaction_gene(reaction_ID,model_ID)''')
            self.database.execute('''CREATE INDEX reactionprotein_ind ON
                                            reaction_protein(reaction_ID,model_ID)''')
            self.database.execute('''CREATE INDEX cluster_ind ON cluster(cluster_num)''')
            self.database.execute('''CREATE INDEX reaction_reversibility_ind ON
                                            reaction_reversibility(reaction_ID)''')
            if self.inchidb is True:
                self.database.execute('''CREATE INDEX original_db_cpdIDs_ind ON
                                                original_db_cpdIDs(ID, inchi_id)''')
        except sqlite3.OperationalError:
            print ('WARNING: Database already exists, reindex indicies')
            self.database.execute('''REINDEX reactioncompound_ind1''')
            self.database.execute('''REINDEX reactioncompound_ind2''')
            self.database.execute('''REINDEX modelreaction_ind1''')
            self.database.execute('''REINDEX modelreaction_ind2''')
            self.database.execute('''REINDEX modelcompound_ind1''')
            self.database.execute('''REINDEX modelcompound_ind2''')
            self.database.execute('''REINDEX model_ind''')
            self.database.execute('''REINDEX reaction_ind''')
            self.database.execute('''REINDEX compound_ind''')
            self.database.execute('''REINDEX reactiongene_ind''')
            self.database.execute('''REINDEX reactionprotein_ind''')
            self.database.execute('''REINDEX cluster_ind''')
            self.database.execute('''REINDEX reaction_reversibility_ind''')
            if self.inchidb is True:
                self.database.execute('''REINDEX original_db_cpdIDs_ind''')
