from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Sets bounds necessary for a specific taraget compound for pulp and runs glpk'

import re
from copy import deepcopy
from timeit import default_timer as timer
import pulp
from ShortestPath import cycle as cy

def verbose_print(verbose, line):
    if verbose:
        print(line)

class IntergerProgram(object):
    """Sets final constraints and solved integer linear program"""
    def __init__(self, db, limit_reactions, limit_cycle, k_paths, cycle, verbose, time_limit, OUTPUT):
        '''initialize class'''
        self.limit_cycle = limit_cycle
        self.limit_reactions = limit_reactions
        self.k_paths = k_paths
        self.cycle = cycle
        self.time_limit = time_limit
        self.verbose = verbose
        self.CYCLE = cy.CycleCheck(db, verbose)
        self.OUTPUT = OUTPUT

    def set_row_bounds(self, lp):
        '''Set row bounds'''
        verbose_print(self.verbose, 'STATUS: Setting compound constraints ...')
        for count, key in enumerate(lp.constraints):
            try:
                inmet = self.incpds.index(self.allcpds[count])
                lp.constraints[key].changeRHS(-10000000000000)
            except ValueError:
                if self.allcpds[count] == self.target:
                    lp.constraints[key].sense = 0
                    lp.constraints[key].changeRHS(1)
                else:
                    lp.constraints[key].changeRHS(0)
            except IndexError:
                pass
        return lp

    def set_objective_function(self, variables):
        '''Set objective function '''
        verbose_print(self.verbose, 'STATUS: Generating objective function coefficients ...')
        start = timer()
        obj = [0]*len(variables)
        for count, variable in enumerate(variables):
            reaction = deepcopy(self.allrxnsrev_dict[str(variable)])
            reaction = re.sub('_F$', '', reaction)
            reaction = re.sub('_R$', '', reaction)
            if reaction not in self.inrxns:
                obj[count] = 1
        end = timer()
        verbose_print(self.verbose, "Time(seconds) to set objective function "+str(self.target)+' '+str(end - start))
        if self.OUTPUT:
            self.OUTPUT.output_timer('Setting objective functions for {}\t{}\t{}\n'.format(self.target, (end-start), (end-start)/60))
        return obj
    
    def set_objective_function_internal(self, variables):
        '''Set objective function '''
        verbose_print(self.verbose, 'STATUS: Generating objective function coefficients (internal)...')
        start = timer()
        obj = [0]*len(variables)
        for count in range(0, len(variables)):
            obj[count] = 1
        end = timer()
        verbose_print(self.verbose, "Time(seconds) to set internal objective function "+str(self.target)+' '+str(end - start))
        if self.OUTPUT:
            self.OUTPUT.output_timer('Setting internal objective function for {}\t{}\t{}\n'.format(self.target, (end-start), (end-start)/60))
        return obj

    def fill_allsolutions(self, solution):
        if solution not in self.allsolutions:
            self.allsolutions.append(solution)

    def ip_calculate(self, lp, variables, obj):
        '''Run glpk solver'''
        start = timer()
        solution = []
        solution_internal = []
        lp.setObjective(pulp.lpSum(obj[i]*variables[i] for i in range(len(obj))))

        verbose_print(self.verbose, 'STATUS: Solving problem for {}...'.format(self.target))
        if self.time_limit == 'None':
            lp.solve(pulp.GLPK(msg=0))
        else:
            tmlim = str(int(self.time_limit)*60)
            lp.solve(pulp.GLPK(msg=0, options=['--tmlim', tmlim]))
        for variable in variables:
            if variable.value() != 0 and variable.value() is not None:
                if variable.name.startswith('Cycle Variable') or variable.name.startswith('cycle'):
                    pass
                else:
                    rxn_sol = deepcopy(self.allrxnsrev_dict[variable.name])
                    if rxn_sol.endswith('_F'):
                        rxn_sol = re.sub('_F$', '', rxn_sol)
                    elif rxn_sol.endswith('_R'):
                        rxn_sol = re.sub('_R$', '', rxn_sol)
                    if rxn_sol not in self.inrxns:
                        solution.append(self.allrxnsrev_dict[variable.name])
                    else:
                        solution_internal.append(self.allrxnsrev_dict[variable.name])
        if self.limit_reactions != 'None' and len(solution) > int(self.limit_reactions):
            print ('STATUS: Path contains too many reaction steps, consider increasing limit up from {} for target {}'.format(self.limit_reactions, self.target))
            solution = []
            solution_internal = []
        end = timer()
        verbose_print(self.verbose, "Time(seconds) to solve for specific path "+str(self.target)+' '+str(end - start))
        if self.OUTPUT:
            self.OUTPUT.output_timer('Solve for specific path for {}\t{}\t{}\n'.format(self.target, (end-start), (end-start)/60))
        return solution, solution_internal

    def run_glpk(self, LP, incpds, inrxns, target_compound_ID, multiplesolutions=True):
        '''Final set up and solve integer linear program'''
        '''Set initial variables'''
        self.inrxns = inrxns
        self.incpds = incpds
        self.target = target_compound_ID
        self.LP = LP
        self.allcpds = self.LP.allcpds
        self.allrxnsrev_dict = self.LP.allrxnsrev_dict
        self.allrxnsrev_dict_rev = self.LP.allrxnsrev_dict_rev
        lp = self.LP.lp
        self.allsolutions = []
        self.k_bounds = []
        self.multiplesolutions = multiplesolutions
        optimalsolutions = []
        optimalsolutions_internal = []
        self.allcyclesolutions = []
        self.solution_threshold = 150
        self.variables_strings = [str(variable) for variable in self.LP.variables]

        '''Set problem bounds and solve'''
        lp = self.set_row_bounds(lp)
        obj = self.set_objective_function(self.LP.variables)
        solution, solution_internal = self.ip_calculate(lp, self.LP.variables, obj)
        self.fill_allsolutions(solution)
        if self.cycle == 'True' and solution:
            '''Check for cycles in pathway'''
            solution, solution_internal, lp, self.LP.variables, obj = self.cycle_constraints(lp, self.LP.variables,
                                                                                             solution, solution_internal,
                                                                                             obj, 0, initialcheck=True)
        if solution_internal and solution:
            self.fill_allsolutions(solution)
            optimalsolutions_internal.append(solution)

        elif not solution_internal and solution:
            self.fill_allsolutions(solution)
            optimalsolutions.append(solution)            

        if multiplesolutions == 'True' and solution:
            '''Check for multiple solutions'''
            optimalsolutions = self.multiple_optimal_solution(lp, self.LP.variables,
                                                              obj, solution, optimalsolutions,
                                                              optimalsolutions_internal, 0)

        if multiplesolutions == 'False' and solution_internal:
            optimalsolutions, variables = self.identify_internal_rxns(self.LP.variables, optimalsolutions,
                                                                      optimalsolutions_internal, lp)

        if len(optimalsolutions) > self.solution_threshold:
            print ('STATUS: Number of solutions {} exceeded limit {} therefore stopping search for target {}'.format(len(optimalsolutions), self.solution_threshold, self.target))
        return optimalsolutions

    def set_weight(self, number_rxn_steps):
        '''retrieve weight to disfavor rxns in path'''
        E = float(1)/float(number_rxn_steps)
        E = round(E, 2)
        E = E-.01
        return E

    def k_number_paths(self, lp, variables, obj, op, op_internal, count_k_paths):
        ''' 'retrieve the next shortest path '''
        if count_k_paths <= self.k_paths:
            obj = self.set_objective_function(variables)
            verbose_print(self.verbose, 'Finding {} suboptimal shortest path'.format(count_k_paths))
            '''Solve for new pathway with new calculated weights'''

            for count_solution, solution in enumerate(self.allsolutions):
                if count_solution not in self.k_bounds:
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append(self.variables_strings.index(self.allrxnsrev_dict_rev[r]))
                    lp += pulp.LpConstraint(pulp.lpSum(1*variables[j] for j in temp),
                                            name='K pathway constraint '+str(count_solution),
                                            sense=-1, rhs=len(temp)-1)
                    self.k_bounds.append(count_solution)
            solution, solution_internal = self.ip_calculate(lp, variables, obj)
            self.fill_allsolutions(solution)

            '''If pathway is greater in steps than the previous path add too optimal solutions
                check for multiple solution'''
            if self.cycle == 'True':
                solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables, solution,
                                                                                         solution_internal, obj, 0)
            if solution_internal:
                if solution and solution not in op_internal+op:
                    op_internal.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj, solution,
                                                        op, op_internal, count_k_paths)
                else:
                    op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
            else:
                if solution and solution not in op_internal+op and len(op) <= self.solution_threshold:
                    op.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj, solution,
                                                        op, op_internal, count_k_paths)
        return op
 
    def identify_internal_rxns(self, variables, op, op_internal, lp):
        ignorerxns = []
        for orig_solution in op_internal:
            for count, variable in enumerate(variables):
                rxn = self.allrxnsrev_dict[str(variable)]
                reaction = deepcopy(self.allrxnsrev_dict[str(variable)])
                reaction = re.sub('_F$', '', reaction)
                reaction = re.sub('_R$', '', reaction)
                if rxn in orig_solution:
                    variable.lowBound = 1
                    variable.upBound = 1
                elif reaction not in self.inrxns:
                    if variable.lowBound == 0 and variable.upBound == 0:
                        ignorerxns.append(rxn) 
                    variable.lowBound = 0
                    variable.upBound = 0

            obj = self.set_objective_function_internal(variables)
            solution, solution_internal = self.ip_calculate(lp, variables, obj)
            if self.cycle == 'True' and solution:
                '''Check for cycles in pathway'''
                solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables, solution+solution_internal,
                                                                               obj, 0, initialcheck=True)
                if solution and len(op) <= self.solution_threshold:
                    op.append(solution)
            else:
                solution = solution + solution_internal
                if solution and len(op) <= self.solution_threshold:
                    op.append(solution)

            '''Add pathway to all optimal solutions'''

            if solution and len(op) <= self.solution_threshold and self.multiplesolutions == 'True':
                verbose_print(self.verbose, 'STATUS: Checking for multiple optimal solutions internal ... ')
                op = self.multiple_optimal_solution_internal(lp, variables, obj, solution, op)

            for variable in variables:
                rxn = self.allrxnsrev_dict[str(variable)]
                if rxn.startswith('EX_') or rxn in ignorerxns:
                    variable.lowBound = 0
                    variable.upBound = 0

                else:
                    variable.lowBound = 0
                    variable.upBound = 1
        return op, variables

    def multiple_optimal_solution_internal(self, lp, variables, obj, originalsolution, op):
        '''Find multiple solutions that use organisms internal reactions'''
        E = self.set_weight(len(originalsolution))
        for rxn in originalsolution:
            index = self.variables_strings.index(self.allrxnsrev_dict_rev[rxn])
            obj[index] = (1+E)
        solution, solution_internal = self.ip_calculate(lp, variables, obj)
        if self.cycle =='True' and solution_internal:
            solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables, 
                                                                           solution+solution_internal,
                                                                           obj, 0)
        else:
            solution = solution+solution_internal
        if len(solution) == len(originalsolution) and len(op) <= self.solution_threshold:
            if solution and solution not in op:
                op.append(solution)
                op = self.multiple_optimal_solution_internal(lp, variables, obj, solution, op)
        return op

    def multiple_optimal_solution(self, lp, variables, obj, originalsolution, op, op_internal, count_k_paths):
        '''Identify multiple solutions'''
        '''Get weight'''
        E = self.set_weight(len(originalsolution))

        '''Set weight for reactions in objective function'''
        for rxn in originalsolution:
            index = self.variables_strings.index(self.allrxnsrev_dict_rev[rxn])
            obj[index] = (1+E)
        solution, solution_internal = self.ip_calculate(lp, variables, obj)
        self.fill_allsolutions(solution)
        '''Check for cycles'''
        if self.cycle == 'True':
            solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables, solution,
                                                                                     solution_internal, obj, 0)

        if len(solution) == len(originalsolution) and len(op) <= self.solution_threshold and len(op_internal) <= self.solution_threshold:
            '''If pathway not already identified add to total solution and check for more '''
            if solution_internal:
                if solution and solution not in op_internal:
                    op_internal.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj, solution,
                                                        op, op_internal, count_k_paths)
                else:
                   op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                   count_k_paths += 1
                   op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
            else:
                if solution and solution not in op_internal+op:
                    op.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj, solution,
                                                        op, [], count_k_paths)

                else:
                    '''If pathway already in solution check for next shortest pathway'''
                    if op_internal:
                        op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                    count_k_paths += 1
                    op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
        else:
            '''If solution is 0 or not the same length as previous solution check for next
                shortest solution'''

            if len(op) <= self.solution_threshold and len(op_internal) <= self.solution_threshold:
                if op_internal:
                    op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                count_k_paths += 1
                op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
            else:
                print ('STATUS: Either solutions types internal {} and/or external {} have exceeded the solution threshold {} for target {}'.format(len(op_internal), len(op), self.solution_threshold, self.target))
                if op_internal:
                    op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
        return op

    def cycle_constraints(self, lp, variables, solution, solution_internal, obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints and resolve
        if cycle is identified
        '''
        verbose_print(self.verbose, 'STATUS: Checking for cycles in the identified pathways for target {}'.format(self.target))

        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ('STATUS: Exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {} for target {}'.format(self.limit_cycle, self.target))
                return ([], [], lp, variables, obj)

            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''
                if solution not in self.allcyclesolutions:
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append(self.variables_strings.index(self.allrxnsrev_dict_rev[r]))
                    self.allcyclesolutions.append(solution)
                    lp += pulp.LpConstraint(pulp.lpSum(1*variables[j] for j in temp), name='cycle pathway constraint '+str(len(self.allcyclesolutions)), sense=-1, rhs=len(temp)-1)
                    solution, solution_internal = self.ip_calculate(lp, variables, obj)

                    if initialcheck is False:
                        if len(solution) == len(original_solution):
                            solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables, solution,
                                                                                                     solution_internal, obj, cycle_count)
                            return (solution, solution_internal, lp, variables, obj)
                        else:
                            print ('WARNING: new solution longer than original {}  for target {}'.format(solution, self.target))
                            return ([], [], lp, variables, obj)
                    else:
                        solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables, solution,
                                                                                                 solution_internal, obj, cycle_count,
                                                                                                 initialcheck=True)
                        return (solution, solution_internal, lp, variables, obj)
                else:
                    return ([], [], lp, variables, obj)

        else:
            return (solution, solution_internal, lp, variables, obj)

    def cycle_constraints_internal(self, lp, variables, solution,  obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints and resolve
        if cycle is identified
        '''
        verbose_print(self.verbose, 'STATUS: Checking for cycles in the identified internal pathways for target {}'.format(self.target))

        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ('STATUS: Path (internal) exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {} for target {}'.format(self.limit_cycle, self.target))
                return ([], lp, variables, obj)

            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''
                if solution not in self.allcyclesolutions:
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append(self.variables_strings.index(self.allrxnsrev_dict_rev[r]))
                    self.allcyclesolutions.append(solution)
                    lp += pulp.LpConstraint(pulp.lpSum(1*variables[j] for j in temp), name='cycle pathway constraint '+str(len(self.allcyclesolutions)), sense=-1, rhs=len(temp)-1)
                    solution, solution_internal = self.ip_calculate(lp, variables, obj)

                    if initialcheck is False:
                        if len(solution+solution_internal) == len(original_solution):
                            solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables,
                                                                                           solution+solution_internal,
                                                                                           obj, cycle_count)
                            return (solution, lp, variables, obj)
                        else:
                            print ('WARNING: new solution longer than original (internal) {} for target {}'.format(solution, self.target))
                            return ([], lp, variables, obj)
                    else:
                        solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables,
                                                                                       solution+solution_internal,
                                                                                       obj, cycle_count, initialcheck=True)
                        return (solution, lp, variables, obj)
                else:
                    return ([], lp, variables, obj) 
        else:
            return (solution, lp, variables, obj)

