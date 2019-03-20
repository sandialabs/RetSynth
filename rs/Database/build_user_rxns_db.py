from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Add user specified reactions to database'

import sqlite3
import re
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

class AddUserRxns2DB(object):
	def __init__(self, DBPath, file_name, model_id='UserAdded', rxntype='bio'):
		self.cnx = sqlite3.connect(DBPath)
		self.file_name = file_name
		self.DB = Q.Connector(DBPath)
		self.DBCPDS = self.DB.get_all_compounds()
		self.model_id = model_id
		self.rxntype = rxntype
		self.newcpds = []
		self.modelcpds = []
		self.newrxns = []
		self.modelrxns = []
		self.rxn_reversibility = []
		self.rxn_cpds = []
		self.rxn_genes = []
		self.rxn_protein = []
		self.IN = indigo.Indigo()
		self.INCHI = indigo_inchi.IndigoInchi(self.IN)
		self.open_user_file()
		self.add_data_2_db()

	def get_rxn_components(self, rxn, ids):
		print (rxn)
		rxn_comp = rxn.split(' -> ')
		rxn_react = rxn_comp[0].split(' + ')
		rxn_prod = rxn_comp[1].split(' + ')
		if ids:
			rxn_react, react_stoic = self.get_stoichometry(rxn_react)
			rxn_prod, prod_stoic = self.get_stoichometry(rxn_prod)
			return(rxn_react, react_stoic, rxn_prod, prod_stoic)
		else:
			return(rxn_react, rxn_prod)

	def get_stoichometry(self, cpds):
		stoichiometry = []
		newcpds = []
		for cpd in cpds:
			#cpd = re.sub('\(|\)', '', cpd)
			match = re.search(r'^\d+\.*\d*', cpd)
			if match:
				print (match.group(0))
				stoichiometry.append(match.group(0))
			else:
				stoichiometry.append(1)
			cpd = re.sub(r'^\d+\.*\d*', '', cpd)
			newcpds.append(cpd)
		return(newcpds, stoichiometry)

	def get_fp_cf(self, cpd):
		if cpd.startswith('InChI'):
			mol = self.INCHI.loadMolecule(cpd)
			cf = mol.grossFormula()
			cf = re.sub(' ', '', cf)
			# fp = mol.fingerprint('full')
			# buffer = fp.toBuffer()
			# buffer_array = [str(i) for i in buffer]
			# buffer_string = ','.join(buffer_array)
			# fp = buffer_string
		else:
			cf = 'None'
			# fp = 'None'
		return(cf)

	def check_cpd_in_db(self,rxn_comps, rxn_stoich, rxn_comps_name, rxn_id, is_prod):
		for count, cpd in enumerate(rxn_comps):
			match = re.search('\_(\w{1})$', cpd)
			if match:
				tmp_cpd = re.sub('\_(\w{1})$', '', cpd)
				cf = self.get_fp_cf(tmp_cpd)
				compartment = match.group(0)
				if cpd not in self.DBCPDS and (cpd, rxn_comps_name[count], compartment, 'None') not in self.newcpds:
					self.newcpds.append((cpd, rxn_comps_name[count], compartment, 'None', cf))
				if (cpd, self.model_id) not in self.modelcpds:
					self.modelcpds.append((cpd, self.model_id))
				if (rxn_id, cpd, 1, 1) not in self.rxn_cpds:
					self.rxn_cpds.append(rxn_id, cpd, is_prod, rxn_stoich[count], 0)
			else:
				compartment = 'c0'
				cf = self.get_fp_cf(cpd)
				if cpd not in self.DBCPDS and (cpd+'_'+'c0', rxn_comps_name[count], compartment, 'None') not in self.newcpds:
					self.newcpds.append((cpd+'_'+'c0', rxn_comps_name[count], compartment, 'None', cf))
				if (cpd+'_'+compartment, self.model_id) not in self.modelcpds:
					self.modelcpds.append((cpd+'_'+compartment, self.model_id))
				if (rxn_id, cpd+'_'+compartment, 1, 1) not in self.rxn_cpds:
					self.rxn_cpds.append((rxn_id, cpd+'_'+compartment, is_prod,  rxn_stoich[count], 0))
	
	def open_user_file(self):
		with open(self.file_name) as fin:
			header = fin.readline()
			for line in fin:
				line = line.strip()
				larray = line.split('\t')
				rxn_react_id, rxn_react_stoic, rxn_prod_id, rxn_prod_stoic = self.get_rxn_components(larray[2], True)
				try:
					if larray[3] == '':
						rxn_react_name = ['None'] * len(rxn_react_id)
						rxn_prod_name = ['None'] * len(rxn_prod_id)
				except IndexError:
					rxn_react_name = ['None'] * len(rxn_react_id)
					rxn_prod_name = ['None'] * len(rxn_prod_id)
				if larray[3] != '':
					rxn_react_name, rxn_prod_name = self.get_rxn_components(larray[3], False)
				else:
					rxn_react_name = ['None'] * len(rxn_react_id)
					rxn_prod_name = ['None'] * len(rxn_prod_id)
				self.check_cpd_in_db(rxn_react_id, rxn_react_stoic, rxn_react_name, larray[0], 0)
				self.check_cpd_in_db(rxn_prod_id, rxn_prod_stoic, rxn_prod_name, larray[0], 1)
				self.newrxns.append((larray[0], larray[1], 'None', self.rxntype))
				try:
					reversibility = str(larray[4]).lower()
				except IndexError:
					reversibility = 'true'
				self.modelrxns.append((larray[0], self.model_id, reversibility))
				self.rxn_reversibility.append((larray[0], reversibility))
				try:
					if larray[5] == '':
						self.rxn_genes.append((larray[0], self.model_id, 'None'))
					else:
						self.rxn_genes.append((larray[0], self.model_id, larray[5]))
				except IndexError:
					self.rxn_genes.append((larray[0], self.model_id, 'None'))
				try:
					if larray[6] == '':
						self.rxn_protein.append((larray[0], self.model_id, 'None'))
					else:
						self.rxn_protein.append((larray[0], self.model_id, larray[6]))
				except IndexError:
					self.rxn_protein.append((larray[0], self.model_id, 'None'))

	def add_data_2_db(self):
		Q = self.cnx.execute('''INSERT INTO model VALUES (?,?)''', (self.model_id, self.file_name))
		Q = self.cnx.execute('''SELECT DISTINCT cluster_num FROM cluster''')
		hits = Q.fetchall()
		uniq_clusters = [i[0] for i in hits]
		self.cnx.execute("INSERT INTO cluster VALUES (?,?)", (len(uniq_clusters)+1, self.model_id))
		self.cnx.commit()
		self.cnx.executemany("INSERT INTO model_compound VALUES (?,?)", self.modelcpds)
		self.cnx.executemany("INSERT INTO model_reaction VALUES (?,?,?)", self.modelrxns)
		self.cnx.executemany("INSERT INTO reaction_gene VALUES (?,?,?)", self.rxn_genes)
		self.cnx.executemany("INSERT INTO reaction_protein VALUES (?,?,?)", self.rxn_protein)
		self.cnx.executemany("INSERT INTO compound VALUES (?, ?, ?, ?, ?)", self.newcpds)
		self.cnx.executemany("INSERT INTO reaction VALUES (?, ?, ?, ?)", self.newrxns)
		self.cnx.executemany("INSERT INTO reaction_compound VALUES (?, ?, ?, ?, ?)", self.rxn_cpds)
		self.cnx.executemany("INSERT INTO reaction_reversibility VALUES (?, ?)", self.rxn_reversibility)
		self.cnx.commit()