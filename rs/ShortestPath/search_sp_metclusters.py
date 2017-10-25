from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Finds shortest path for all unique metabolic clusters'

class SearchMetabolicClusters(object):
    """
    Species with the same metabolism are grouped into a defined cluster in the database,
    this goes through all defined clusters and determines the number of
    reactions that need to be added to each cluster to produce a target compound
    """

    def __init__(self, target_compound_ID, LP, IP, OUTPUT, db):
        '''Intialize class'''
        self.OUTPUT = OUTPUT
        self.DB = db
        self.LP = LP
        self.IP = IP
        self.target = target_compound_ID
        cpdname = self.DB.get_compound_name(self.target)
        self.OUTPUT.output_open_paths_all_organism_file(cpdname)
        self.get_sp_4_clusters()

    def get_sp_4_clusters(self):
        '''
        Goes through each metabolic cluster and identifies number of reactions that need to be
        added to produce a target compound
        '''
        self.total_sp_clusters = []
        uniq_clusters = self.DB.get_uniq_metabolic_clusters()
        for count, cluster_num in enumerate(uniq_clusters):
            orgs = self.DB.get_models_from_cluster(cluster_num)
            org_names = []
            for o in orgs:
                name = self.DB.get_organism_name(o)
                org_names.append(name)
            model = orgs[0]
            inmets = self.DB.get_compounds_in_model(model)
            inrxns = self.DB.get_reactions_in_model(model)
            optimal_pathways = self.IP.run_glpk(self.LP, inmets, inrxns, self.target,
                                                multiplesolutions=False)
            self.total_sp_clusters.append([self.target, len(optimal_pathways[0]), orgs])
            self.OUTPUT.output_paths_all_organisms(self.target, len(optimal_pathways[0]),
                                                   orgs, org_names)
            print ('{} out of {} unique metabolic clusters shortest paths identified'.format(count+1, len(uniq_clusters)))
