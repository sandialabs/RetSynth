from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Translates kbase genes in database'

import re
import os
PATH = os.path.dirname(os.path.abspath(__file__))

class TranslateKbaseGenes(object):
    """Translate kbase gene IDs to there names"""
    def __init__(self, DB, file_gene_translate, OUTPUTPATH=PATH):
        self.DB = DB
        self.OUTPUTPATH = OUTPUTPATH
        self.file_gene_translate = PATH + '/' + file_gene_translate
        self.translation_dict = {}
        if os.path.isfile(OUTPUTPATH + '/GeneTranslations.txt') is False:
            print ('STATUS: Generating new gene translation file ...')
            self.translate_genes()
            self.open_translation_file()

        else:
            print ('STATUS: Loading gene translations from preexisting file')
            gene_translatefile = open(OUTPUTPATH + '/GeneTranslations.txt')
            for line in gene_translatefile:
                larray = line.strip().split()
                self.translation_dict[larray[0]] = larray[1]

    def open_translation_file(self):
        '''
        Opens file that has translation information and processes the file
        '''
        with open(self.OUTPUTPATH + '/GeneTranslations.txt', 'w') as outputfile:
            print ('STATUS: Retrieving translations of genes ...')
            with open(self.file_gene_translate) as file_name:
                kbase_genes = {}
                for count, line in enumerate(file_name):
                    larray = line.strip('\n').split('   ')
                    kbase_genes[larray[0]] = larray[1]
                intersect = set(kbase_genes.keys()).intersection(self.DBgenes)
                for i in intersect:
                    self.translation_dict[i] = kbase_genes[i]
                    outputfile.write(i + '\t' + kbase_genes[i] + '\n')
                if len(self.translation_dict) != len(self.DBgenes):
                    print ('WARNING: discrepency between genes in DB and genes retrieved from translation file')
            print ('STATUS: Finished retrieving translations ...')

    def translate_genes(self):
        '''
        Retrieves all the genes in the database whose name needs to be
        identified
        '''
        print ('STATUS: Retrieving all genes in database ...')
        models = self.DB.get_all_models()
        self.DBgenes = []
        for count, model in enumerate(models):
            print ('{} out of {} models processed'.format(count+1, len(models)))
            rxns = self.DB.get_reactions_in_model(model[0])
            for rxn in rxns:
                genes = self.DB.get_genes(rxn, model[0])
                genes = re.sub('\)', '', genes)
                genes = re.sub('\(', '', genes)
                individualgenes = genes.strip().split()
                for g in individualgenes:
                    if g == 'and' or g == 'or' or g == 'None' or g == 'unknown' or g == 'Unknown':
                        pass
                    else:
                        self.DBgenes.append(g)
        self.DBgenes = list(set(self.DBgenes))
        print ('STATUS: Finished retrieving all genes in database')
