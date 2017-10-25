from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Generates Sqlite database'

import sqlite3
from Database import build_tablesMP as bt

class Createdb(object):
    """Generates new sqlite database and creates table and indicies"""
    def __init__(self, db, dumpfile_dir, inchidb, rxntype):
        self.df_dir = dumpfile_dir
        self.inchidb = inchidb
        self.sqlite_database = sqlite3.connect(db)
        self.db = db
        self.rxntype = rxntype
        self.create_model_tbls()

    def create_model_tbls(self):
        '''
        Creates tables and indicies in SQLite Database
        '''
        try:
            self.sqlite_database.execute('''CREATE table model
                                            (ID text,file_name text)''')
            self.sqlite_database.execute('''CREATE table compound
                                            (ID text, name text,compartment text)''')
            self.sqlite_database.execute('''CREATE table compartments (ID text, name text)''')
            self.sqlite_database.execute('''CREATE table model_compound
                                            (cpd_ID text, model_ID text)''')
            self.sqlite_database.execute('''CREATE table reaction
                                            (ID text, name text, type text)''')
            self.sqlite_database.execute('''CREATE table model_reaction
                                            (reaction_ID text, model_ID text,is_rev bit(1))''')
            self.sqlite_database.execute('''CREATE table reaction_compound
                                            (reaction_ID text, cpd_ID text, is_prod bit(1), 
                                            stoichiometry int, filenum int)''')
            self.sqlite_database.execute('''CREATE table reaction_reversibility
                                            (reaction_ID text, is_reversible bit(1))''')
            self.sqlite_database.execute('''CREATE table reaction_gene
                                            (reaction_ID text, model_ID text, gene_ID text)''')
            self.sqlite_database.execute('''CREATE table reaction_protein
                                            (reaction_ID text, model_ID text, protein_ID text)''')
            self.sqlite_database.execute('''CREATE table cluster
                                            (cluster_num text, ID text)''')
            if self.inchidb is True:
                self.sqlite_database.execute('''CREATE table original_db_cpdIDs
                                                (ID text, inchi_id text)''')
            print ('STATUS: Generating new database ...')
        except sqlite3.OperationalError:
            print ('WARNING: Database already exists, adding new xml files to existing database ...')
        bt.BuildTables(self.df_dir, self.inchidb, self.db, self.rxntype)
        try:
            self.sqlite_database.execute('''CREATE INDEX reactioncompound_ind1 ON
                                            reaction_compound(reaction_ID,cpd_ID,is_prod)''')
            self.sqlite_database.execute('''CREATE INDEX reactioncompound_ind2 ON
                                            reaction_compound(cpd_ID,is_prod)''')
            self.sqlite_database.execute('''CREATE INDEX modelreaction_ind1 ON
                                            model_reaction(model_ID)''')
            self.sqlite_database.execute('''CREATE INDEX modelreaction_ind2 ON
                                            model_reaction(reaction_ID)''')
            self.sqlite_database.execute('''CREATE INDEX modelcompound_ind1 ON
                                            model_compound(model_ID)''')
            self.sqlite_database.execute('''CREATE INDEX modelcompound_ind2 ON
                                            model_compound(cpd_ID)''')
            self.sqlite_database.execute('''CREATE INDEX model_ind ON
                                            model(ID)''')
            self.sqlite_database.execute('''CREATE INDEX reaction_ind ON
                                            reaction(ID)''')
            self.sqlite_database.execute('''CREATE INDEX compound_ind ON
                                            compound(ID)''')
            self.sqlite_database.execute('''CREATE INDEX reactiongene_ind ON
                                            reaction_gene(reaction_ID,model_ID)''')
            self.sqlite_database.execute('''CREATE INDEX reactionprotein_ind ON
                                            reaction_protein(reaction_ID,model_ID)''')
            self.sqlite_database.execute('''CREATE INDEX cluster_ind ON cluster(cluster_num)''')
            self.sqlite_database.execute('''CREATE INDEX reaction_reversibility_ind ON
                                            reaction_reversibility(reaction_ID)''')
            if self.inchidb is True:
                self.sqlite_database.execute('''CREATE INDEX original_db_cpdIDs_ind ON
                                                original_db_cpdIDs(ID, inchi_id)''')
        except sqlite3.OperationalError:
            print ('WARNING: Database already exists, reindex indicies')
            self.sqlite_database.execute('''REINDEX reactioncompound_ind1''')
            self.sqlite_database.execute('''REINDEX reactioncompound_ind2''')
            self.sqlite_database.execute('''REINDEX modelreaction_ind1''')
            self.sqlite_database.execute('''REINDEX modelreaction_ind2''')
            self.sqlite_database.execute('''REINDEX modelcompound_ind1''')
            self.sqlite_database.execute('''REINDEX modelcompound_ind2''')
            self.sqlite_database.execute('''REINDEX model_ind''')
            self.sqlite_database.execute('''REINDEX reaction_ind''')
            self.sqlite_database.execute('''REINDEX compound_ind''')
            self.sqlite_database.execute('''REINDEX reactiongene_ind''')
            self.sqlite_database.execute('''REINDEX reactionprotein_ind''')
            self.sqlite_database.execute('''REINDEX cluster_ind''')
            self.sqlite_database.execute('''REINDEX reaction_reversibility_ind''')
            if self.inchidb is True:
                self.sqlite_database.execute('''REINDEX original_db_cpdIDs_ind''')
