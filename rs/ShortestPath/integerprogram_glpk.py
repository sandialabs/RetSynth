from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Sets bounds necessary for a specific taraget compound for pyglpk and runs glpk'

from copy import deepcopy
import re
from ShortestPath import cycle as cy
from tqdm import tqdm

class IntergerProgram(object):
    """Sets final constraints and solved integer linear program"""
    def __init__(self, db, limit_reactions, limit_cycle, k_paths, cycle):
        '''initalize class'''
        self.limit_cycle = limit_cycle
        self.limit_reactions = limit_reactions
        self.k_paths = k_paths
        self.cycle = cycle
        self.CYCLE = cy.CycleCheck(db)
        self.weight_dict = {}

    def set_row_bounds(self, lp):
        '''Set row bounds'''
        print ('STATUS: Setting compound constraints ...')
        for count, r in enumerate(lp.rows):
            try:
                inmet = self.incpds.index(self.allcpds[count])
                r.bounds = None, None
            except ValueError:
                met = self.allcpds[count]
                if met == self.target:
                    r.bounds = 1
                else:
                    r.bounds = 0, None
            except IndexError:
                pass
        return lp

    def set_objective_function(self, lp):
        '''Set objective function '''
        print ('STATUS: Generating objective function coefficients ...')
        obj = [0]*len(lp.cols)
        for count, r in enumerate(self.allrxnsrev):
            reaction = deepcopy(r)
            reaction = re.sub('_F$', '', reaction)
            reaction = re.sub('_R$', '', reaction)
            if reaction not in self.inrxns:
                obj[count] = 1
        return obj

    def set_objective_function_internal(self, lp):
        '''Set objective function '''
        print ('STATUS: Generating objective function coefficients (internal)...')
        obj = [1]*len(lp.cols)
        for count, r in enumerate(self.allrxnsrev):
            obj[count] = 1

        return obj

    def fill_allsolutions(self, solution):
        if solution not in self.allsolutions:
            self.allsolutions.append(solution)

    def ip_calculate(self, lp, obj):
        '''Run glpk solver'''
        solution = []
        solution_internal = []
        lp.obj[:] = obj
        lp.intopt()
        if lp.status == 'undef':
            '''Check if solution was found'''
            pass
        else:
            if lp.obj.value > 0:
                for c in lp.cols:
                    if c.primal != 0:
                        if c.name == 'Cycle Variable':
                            pass
                        else:
                            rxn_sol = self.allrxnsrev[self.LP.rxnnames.index(c.name)]
                            if rxn_sol.endswith('_F'):
                                rxn_sol = re.sub('_F$', '', rxn_sol)
                            elif rxn_sol.endswith('_R'):
                                rxn_sol = re.sub('_R$', '', rxn_sol)
                            if rxn_sol not in self.inrxns:
                                solution.append(self.allrxnsrev[self.LP.rxnnames.index(c.name)])
                            else:
                                solution_internal.append(self.allrxnsrev[self.LP.rxnnames.index(c.name)])
                if self.limit_reactions != 'None' and len(solution) > int(self.limit_reactions):
                    print ('STATUS: Path contains too many reaction steps, consider increasing limit up from {} first line'.format(self.limit_reactions))
                    solution = []
                    solution_internal = []
            else:
                solution = []
                solution_internal = []
        return solution, solution_internal

    def run_glpk(self, LP, incpds, inrxns, target_compound_ID,
                 multiplesolutions, eliminatedrxns=False):
        '''Final set up and solve integer linear program'''

        '''Set initial variables'''
        self.inrxns = inrxns
        self.incpds = incpds
        self.target = target_compound_ID
        self.LP = LP
        self.allrxnsrev = self.LP.allrxnsrev
        self.allcpds = self.LP.allcpds
        self.k_bounds = []
        self.allsolutions = []
        lp = self.LP.lp
        optimalsolutions = []
        optimalsolutions_internal = []
        self.allcyclesolutions = []

        '''Set problem bounds and solve'''
        lp = self.set_row_bounds(lp)
        obj = self.set_objective_function(lp)
        solution, solution_internal = self.ip_calculate(lp, obj)
        self.fill_allsolutions(solution)

        if self.cycle == 'True' and solution:
            '''Check for cycles in pathway'''
            solution, solution_internal, lp, obj = self.cycle_constraints(lp, solution, solution_internal,
                                                                          obj, 0, initialcheck=True)

        '''Add pathway to all optimal solutions'''
        if solution_internal:
            self.fill_allsolutions(solution)
            if solution:
                optimalsolutions_internal.append(solution)
        else:
            self.fill_allsolutions(solution)
            if solution:
                optimalsolutions.append(solution)

        if multiplesolutions == 'True' and solution:
            '''Check for multiple solutions'''
            print ('STATUS: Checking for multiple optimal solutions ... ')
            optimalsolutions = self.multiple_optimal_solution(lp, obj, solution,
                                                              optimalsolutions, optimalsolutions_internal, 0)

        return optimalsolutions

    def set_weight(self, number_rxn_steps):
        '''retrieve weight to disfavor rxns in path'''
        E = float(1)/float(number_rxn_steps)
        E = round(E, 2)
        E = E-.01
        return E


    def k_number_paths(self, lp, obj, op, op_internal, count_k_paths):
        ''' 'retrieve the next shortest path '''
        if count_k_paths <= self.k_paths:
            obj = self.set_objective_function(lp)
            print ('Finding {} shortest path'.format(count_k_paths))
            for count_solutions, solution in enumerate(self.allsolutions):
                if count_solutions not in self.k_bounds:
                    lp.rows.add(1)
                    new_row = lp.rows[-1]
                    new_row.name = 'K pathway constraint '+str(count_solutions)
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append((self.allrxnsrev.index(r), 1))
                    new_row.matrix = temp
                    new_row.bounds = 0, len(temp)-1
                    self.k_bounds.append(count_solutions)
            solution, solution_internal = self.ip_calculate(lp, obj)
            self.fill_allsolutions(solution)
            

            '''If pathway is greater in steps than the previous path,
             check for cycles and then add too optimal solutions
             check for multiple solution'''
            if self.cycle == 'True':
                solution, solution_internal, lp, obj = self.cycle_constraints(lp, solution, solution_internal, obj, 0)

            if solution_internal:
                if solution and solution not in op_internal+op:
                    op_internal.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, obj, solution, op, op_internal, count_k_paths)
                else:
                    op, lp = self.identify_internal_rxns(op, op_internal, lp)
            else:
                if solution and solution not in op_internal+op:
                    op.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, obj, solution, op, op_internal, count_k_paths)
 
        return op

    def multiple_optimal_solution(self, lp, obj, originalsolution, op, op_internal, count_k_paths):
        '''Identify multiple solutions'''
        '''Get weight'''
        E = self.set_weight(len(originalsolution))
       
        '''Set weight for reactions in objective function'''
        for rxn in originalsolution:
            index = self.allrxnsrev.index(rxn)
            obj[index] = (1+E)
        solution, solution_internal = self.ip_calculate(lp, obj)
        self.fill_allsolutions(solution)
        '''Check for cycles'''
        if self.cycle == 'True':
            solution, solution_internal, lp, obj = self.cycle_constraints(lp, solution, solution_internal, obj, 0)

        if len(solution) == len(originalsolution):
            ''' Check if solution is of the same number of reaction steps and does
                not have 0 reaction steps'''
            if solution_internal:
                if solution and solution not in op_internal+op:
                    op_internal.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, obj, solution, op, op_internal, count_k_paths)
                else:
                    op, lp = self.identify_internal_rxns(op, op_internal, lp)
                    count_k_paths += 1
                    op = self.k_number_paths(lp, obj, op, [], count_k_paths)
            else:
                if solution and solution not in op_internal+op:
                    '''If pathway not already identified add to total solution and check for more '''
                    op.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, obj, solution, op, op_internal, count_k_paths)

                else:
                        
                    '''If pathway already in solution check for next shortest pathway'''
                    if op_internal:
                        op, lp = self.identify_internal_rxns(op, op_internal, lp)
                    count_k_paths += 1
                    op = self.k_number_paths(lp, obj, op, [], count_k_paths)

        else:
            '''If solution is 0 or not the same length as previous solution check for next
                shortest solution'''
            if op_internal:
                op, lp = self.identify_internal_rxns(op, op_internal, lp)            
            count_k_paths += 1
            op = self.k_number_paths(lp, obj, op, [], count_k_paths)
            print ('STATUS: No more multiple solutions...')
        return op

    def identify_internal_rxns(self, op, op_internal, lp):
        '''Identify the internal reactions necessary to get target compound'''
        print ('STATUS: Identifying internal reactions')
        ignorerxns = []
        for orig_solution in op_internal:
            for count, rxn in enumerate(self.allrxnsrev):
                reaction = deepcopy(rxn)
                reaction = re.sub('_F$', '', reaction)
                reaction = re.sub('_R$', '', reaction)
                if rxn in orig_solution:
                    lp.cols[count].bounds = 1,1
                elif reaction not in self.inrxns:
                    if lp.cols[count].bounds[0] == 0 and lp.cols[count].bounds[1] == 0:
                        ignorerxns.append(rxn) 
                    lp.cols[count].bounds = 0,0

            obj = self.set_objective_function_internal(lp)
            solution, solution_internal = self.ip_calculate(lp, obj)
            if self.cycle == 'True' and solution:
                '''Check for cycles in pathway'''
                solution, lp, obj = self.cycle_constraints_internal(lp, solution+solution_internal, obj, 0, initialcheck=True)
                if solution:
                    op.append(solution)

            else:
                solution = solution + solution_internal
                if solution:
                    op.append(solution)
 
            '''Add pathway to all optimal solutions'''
            if solution:
                print ('STATUS: Checking for multiple optimal solutions internal ... ')
                op = self.multiple_optimal_solution_internal(lp, obj, solution, op)

            for count, rxn in enumerate(self.allrxnsrev):
                if rxn.startswith('EX_') or rxn in ignorerxns:
                    lp.cols[count].bounds = 0, 0
                else:
                    lp.cols[count].bounds = 0, 1
        return op, lp

    def multiple_optimal_solution_internal(self, lp, obj, originalsolution, op):
        E = self.set_weight(len(originalsolution))
        for rxn in originalsolution:
            index = self.allrxnsrev.index(rxn)
            obj[index] = (1+E)
        solution, solution_internal = self.ip_calculate(lp, obj)
        if self.cycle =='True' and solution_internal:
            solution, lp, obj = self.cycle_constraints_internal(lp, solution+solution_internal,
                                                                obj, 0, initialcheck=False)
        else:
            solution = solution+solution_internal

        if len(solution) == len(originalsolution):
            if solution and solution not in op:
                op.append(solution)
                op = self.multiple_optimal_solution_internal(lp, obj, solution, op)
        return op

    def cycle_constraints(self, lp, solution, solution_internal, obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints
        and resolve if cycle is identified
        '''
        print ('STATUS: Checking for cycles in the identified pathway')
        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)
        originalsolution_internal = deepcopy(solution_internal)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ("STATUS: Exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {} ".format(self.limit_cycle))
                return ([], [], lp, obj)
            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''

                if solution not in self.allcyclesolutions:
                    lp.rows.add(1)
                    new_row = lp.rows[-1]
                    new_row.name = 'cycle pathway constraint '+str(len(self.allcyclesolutions))
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append((self.allrxnsrev.index(r), 1))
                    new_row.matrix = temp
                    new_row.bounds = 0, len(temp)-1
                    self.allcyclesolutions.append(solution)
                solution, solution_internal = self.ip_calculate(lp, obj)
                if initialcheck is False:
                    if len(solution) == len(original_solution):
                        solution, solution_internal, lp, obj = self.cycle_constraints(lp, solution, 
                                                                                      solution_internal, obj,
                                                                                      cycle_count)
                        return (solution, solution_internal, lp, obj)

                    else:
                        print ('WARNING: new solution longer than original {}'.format(solution))
                        return ([], [], lp, obj)
                else:
                    solution, solution_internal, lp, obj = self.cycle_constraints(lp, solution, solution_internal,
                                                                        obj, cycle_count,
                                                                        initialcheck=True)
                    return (solution, solution_internal, lp, obj)
        else:
            return (solution, solution_internal, lp, obj)

    def cycle_constraints_internal(self, lp, solution, obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints
        and resolve if cycle is identified for solutions with internal reactions
        '''
        print ('STATUS: Checking for cycles in the identified pathway (internal)')
        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ("STATUS: Exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {} (internal)".format(self.limit_cycle))
                return ([], lp, obj)
            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''
                if solution not in self.allcyclesolutions:
                    lp.rows.add(1)
                    new_row = lp.rows[-1]
                    new_row.name = 'cycle pathway constraint '+str(len(self.allcyclesolutions))
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append((self.allrxnsrev.index(r), 1))
                    new_row.matrix = temp
                    new_row.bounds = 0, len(temp)-1
                    self.allcyclesolutions.append(solution)
                solution, solution_internal = self.ip_calculate(lp, obj)
                if initialcheck is False:
                    if len(solution+solution_internal) == len(original_solution):
                        solution, lp, obj = self.cycle_constraints_internal(lp, solution+solution_internal, obj, cycle_count)
                        return (solution, lp, obj)

                    else:
                        print ('WARNING: new solution longer than original {} (internal)'.format(solution))
                        return ([], lp, obj)
                else:
                    solution, lp, obj = self.cycle_constraints_internal(lp, solution+solution_internal,
                                                                        obj, cycle_count, initialcheck=True)
                    return (solution, lp, obj)
        else:
            return (solution, lp, obj)
  