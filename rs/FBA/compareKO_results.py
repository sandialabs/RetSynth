from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Compare flux balance analysis \
                   results from wildtype to results with added reactions'

class CompareKO(object):
    """
    Compares flux values between simulation with added reactions simulations
    for each reaction knockout
    """
    def __init__(self, target_compound_ID, compounds_dict, wtresults, exresults, ex_paths, DB, fold_threshold=1.5):
        '''Initalize classes'''
        self.wtresults = wtresults
        self.exresults = exresults
        self.target = target_compound_ID
        self.ex_paths = ex_paths
        self.compounds_dict = compounds_dict
        self.DB = DB
        self.fold_threshold = float(fold_threshold)+float(1)
        self.fluxchange = {}
        self.externalrxnfluxes = {}
        self.ex_rxns = []
        for key, path in self.ex_paths.iteritems():
            for r in path:
                if r not in self.ex_rxns:
                    self.ex_rxns.append(r)
        self.get_flux_differences()
        self.external_pathanalyze_fluxes()

    def analyze_fluxes(self, r, rko):
        '''
        Determines if flux differences between simulations is greater than 2.5 fold
        '''
        name = self.DB.get_reaction_name(r)
        if r in self.wtresults.x_dict:
            if r in self.ex_rxns:
                if self.wtresults.x_dict[r] != 0 and self.exresults[rko].x_dict[r] != 0:
                    if self.exresults[rko].x_dict[r] >= self.wtresults.x_dict[r]:
                        self.fluxchange[rko][r] = '\t'.join([name, str(self.wtresults.x_dict[r]),
                                                             str(self.exresults[rko].x_dict[r])])
                elif self.wtresults.x_dict[r] == 0 and self.exresults[rko].x_dict[r] != 0:
                    if self.exresults[rko].x_dict[r] <= -1 or self.exresults[rko].x_dict[r] >= 1:
                        self.fluxchange[rko][r] = '\t'.join([name, str(self.wtresults.x_dict[r]),
                                                             str(self.exresults[rko].x_dict[r])])
            else:
                if self.wtresults.x_dict[r] != 0 and self.exresults[rko].x_dict[r] != 0:
                    upperfold = self.wtresults.x_dict[r]*self.fold_threshold
                    changefold = upperfold-self.wtresults.x_dict[r]
                    changefold = abs(changefold)
                    if self.wtresults.x_dict[r] < 0:
                        lowerfold = self.wtresults.x_dict[r]+changefold
                    elif self.wtresults.x_dict[r] > 0:
                        lowerfold = self.wtresults.x_dict[r]-changefold
                    if upperfold < 0:
                        if (self.exresults[rko].x_dict[r] <= upperfold or
                                self.exresults[rko].x_dict[r] >= lowerfold):
                            self.fluxchange[rko][r] = '\t'.join([name,
                                                                 str(self.wtresults.x_dict[r]),
                                                                 str(self.exresults[rko].x_dict[r])
                                                                ])
                    else:
                        if (self.exresults[rko].x_dict[r] <= lowerfold or
                                self.exresults[rko].x_dict[r] >= upperfold):
                            self.fluxchange[rko][r] = '\t'.join([name,
                                                                 str(self.wtresults.x_dict[r]),
                                                                 str(self.exresults[rko].x_dict[r])
                                                                ])
                elif self.wtresults.x_dict[r] == 0 and self.exresults[rko].x_dict[r] != 0:
                    if self.exresults[rko].x_dict[r] <= -1 or self.exresults[rko].x_dict[r] >= 1:
                        self.fluxchange[rko][r] = '\t'.join([str(self.wtresults.x_dict[r]),
                                                             str(self.exresults[rko].x_dict[r])])

    def get_flux_differences(self):
        '''
        Gets reactions for the flux difference is greater than 2.5 fold
        '''
        for rko in self.exresults.keys():
            self.fluxchange[rko] = {}
            for r in self.exresults[rko].x_dict.keys():
                self.analyze_fluxes(r, rko)

    def external_pathanalyze_fluxes(self):
        '''
        Extracts external path carrying the most flux
        '''
        pathwayflux = {}
        self.objective_function_ko = {}
        self.maxflux = {}
        self.maxpath = {}
        for rko in self.fluxchange:
            pathwayflux[rko] = {}
            self.maxflux[rko] = {}
            self.maxpath[rko] = {}
            if len(self.ex_paths) != 0:
                for key, path in self.ex_paths.iteritems():
                    flux = 0
                    count = 0
                    for r in path:
                        if r in self.fluxchange[rko]:
                            values = self.fluxchange[rko][r].split('\t')
                            flux += abs(float(values[2]))
                        else:
                            count += 1
                    if count == len(path):
                        pathwayflux[rko][key] = 0
                    else:
                        pathwayflux[rko][key] = flux
                self.maxflux[rko] = max(pathwayflux[rko].values())
                self.maxpath[rko] = pathwayflux[rko].keys()[pathwayflux[rko].values().index(self.maxflux[rko])]
            else:
                self.maxpath[rko] = 'No added path'
                self.maxflux[rko] = self.exresults[rko].x_dict['Sink_'+self.target]
            glucose = True
            try:
                glucoseimport = self.exresults[rko].x_dict['EX_cpd00027_e0']
            except KeyError:
                glucose = False
            if glucose:
                objectivesol = self.exresults[rko].x_dict['Sink_'+self.compounds_dict[self.target]]
                glucoseimport = round(glucoseimport, 2)
                try:
                    ty = abs(round(round(objectivesol, 2)/round(glucoseimport,2), 2))
                    self.objective_function_ko[rko] = ty
                except ZeroDivisionError:
                    self.objective_function_ko[rko] = float(0)
            else:
                self.objective_function_ko[rko] = float(0)

