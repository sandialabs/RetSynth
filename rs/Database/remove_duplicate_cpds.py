from __future__ import print_function

__author__ = 'Leanne Whitmore and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov and cmhudso@sandia.gov'
__description__ = 'combined overlapping IDs'
import sqlite3
from tqdm import tqdm
from Database import query as Q

class OverlappingCpdIDs(object):
    def identify_overlapping_ids(self):
        '''find cpd IDs that are have different IDs (i.e. inchi and cpd for the same compound)
            and fix'''
        self.inchicpds2keggids = {}
        self.olids = {}
        self.keggidol = {}
        self.cpds = self.DB.get_all_compounds()
        
        print ('STATUS: get inchi compounds')
        for cpd in self.cpds:
            if cpd.startswith('InChI'):
                self.inchicpds2keggids[cpd] = {}
                self.inchicpds2keggids[cpd]['compartment'] = self.DB.get_compound_compartment(cpd)
                self.inchicpds2keggids[cpd]['keggid'] = self.DB.get_kegg_cpd_ID(cpd)
        print ('STATUS: get overlapping compounds')
        query_search = 'InChI'+str('%')
        query = "select ID from compound where ID not like '%s'" % query_search
        conn, cnx = self.DB.connect_to_database()
        Q = cnx.execute(query)
        hits = Q.fetchall()
        self.keggidsonly = [i[0] for i in hits]

        Q = cnx.execute("SELECT * FROM compartments")
        hits = Q.fetchall()
        compartments =  [i[0] for i in hits]
        for cpd in tqdm(self.inchicpds2keggids):
            if self.inchicpds2keggids[cpd]['keggid']+'_'+self.inchicpds2keggids[cpd]['compartment'] in self.keggidsonly:
                self.keggidol[self.inchicpds2keggids[cpd]['keggid']+'_'+self.inchicpds2keggids[cpd]['compartment']] = cpd

        print ('STATUS: removing duplicate kegg ids from compound table')
        for cpdkeggid in tqdm(self.keggidol):
            cnx.execute("DELETE FROM compound where ID = ?", (cpdkeggid,))
            cnx.execute("UPDATE model_compound SET cpd_ID=? WHERE cpd_ID=?", (self.keggidol[cpdkeggid], cpdkeggid))
            cnx.execute("UPDATE reaction_compound SET cpd_ID=? WHERE cpd_ID=?", (self.keggidol[cpdkeggid], cpdkeggid))
        conn.commit()
    
    def __init__(self, database):
        '''Initialize'''
        self.DB = Q.Connector(database)
        self.identify_overlapping_ids()
