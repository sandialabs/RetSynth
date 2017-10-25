__author__='Leanne Whitmore'
__email__='lwhitmo@sandia.gov'
__description__='Queries compound database'

import sqlite3

class Connector(object):
	"""Connects to a generated database"""
	def __init__(self,db):
		self.cnx=sqlite3.connect(db)

	def get_inchi(self,compound_name):
		'''Retrieve inchis from pubchemDB'''
		if compound_name.strip() == "":
			return str(None)
		query1 = "select ID from synonym where synonym_value like ?"
		Q1=self.cnx.execute(query1, (compound_name,))
		if Q1.arraysize > 0:
			result = Q1.fetchone()
		if result is not None:
			query2 = "select inchi_value from inchi where ID = ?"
			Q2=self.cnx.execute(query2, (result[0],))
			if Q2.arraysize > 0:
				result = Q2.fetchone()
				if result is None:
					return str(None)
				else:
					return(str(result[0])) 
			else:
				return str(None)
		else:
			return str(None)
