from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Sets bounds necessary for a specific taraget compound for pulp and runs glpk'

import re
from copy import deepcopy
from timeit import default_timer as timer
import pulp

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
        self.DB = db
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
        #print (lp.objective.value())
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

    def set_lp_problem(self, lp, variables):
        lp = self.set_row_bounds(lp)
        obj = self.set_objective_function(variables)
        return (lp, obj)

    def initiate_cycle_check(self, solution, solution_internal, lp, variables, obj, cycle_check_count, initialcheck_value=False):
        if self.cycle == 'True' and solution:
            '''Check for cycles in pathway'''
            solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables, solution, solution_internal,
                                                                                     obj, cycle_check_count, initialcheck_value)
        elif not solution:
            verbose_print(self.verbose, 'STATUS: Solution for compound {} is empty therefore not performing cycle check'.format(self.target))
        return (solution, solution_internal, lp, variables, obj)

    def initiate_internal_cycle_check(self, solution, lp, variables, obj, cycle_check_count, length_external):
        if self.cycle == 'True' and solution:
            solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables, solution, obj, cycle_check_count, length_external)
        elif not solution:
            verbose_print(self.verbose, 'STATUS: Solution for compound {} is empty therefore not performing cycle check'.format(self.target))
        return (solution, lp, variables, obj)

    def filling_optimal_solution_arrays(self, solution, solution_internal, op, op_internal):
        check_new_solution=True
        check_solution_threshold=True
        if solution_internal and solution and solution not in op_internal and len(op_internal) <= self.solution_threshold and len(op) <= self.solution_threshold:
            op_internal.append(solution)

        elif not solution_internal and solution and solution not in op and len(op_internal) <= self.solution_threshold and len(op) <= self.solution_threshold:
            op.append(solution)
        
        elif not solution_internal and not solution:
            verbose_print(self.verbose, 'STATUS: Both solution and solution internal are empty for target {} therefore not filling final solution arrays'.format(self.target))
            check_new_solution=False

        elif solution in op or solution in op_internal:
            verbose_print(self.verbose, 'STATUS: Solution for target {} already idenfied therefore not filling final solution arrays'.format(self.target))
            check_new_solution=False

        if len(op_internal) > self.solution_threshold and len(op) > self.solution_threshold:
            check_solution_threshold=False
            check_new_solution=False

        if len(op) > self.solution_threshold:
            check_solution_threshold=False
            check_new_solution=False
 
        return(op, op_internal, check_new_solution, check_solution_threshold)

    def initiate_multiple_solutions(self, solution, lp, variables, obj, op, op_internal, count_k_paths):
        if self.multiplesolutions == 'True' and solution:
            '''Check for multiple solutions'''
            op = self.multiple_optimal_solution(lp, variables, obj, solution, op, op_internal, count_k_paths)
        elif not solution: 
            verbose_print(self.verbose, 'STATUS: Solution for compound {} is empty therefore not performing going to get multiple optimal solutions'.format(self.target))

        return op
    
    # def initiate_identification_of_internal_rxns(self,):

    def run_glpk(self, LP, incpds, inrxns, target_compound_ID, multiplesolutions=True):
        '''Final set up and solve integer linear program'''
        '''Set initial variables'''
        self.inrxns = inrxns
        self.incpds = incpds
        self.target = target_compound_ID
        self.allcpds = LP.allcpds
        self.allrxnsrev_dict = LP.allrxnsrev_dict
        self.allrxnsrev_dict_rev = LP.allrxnsrev_dict_rev
        self.allsolutions = []
        self.k_bounds = []
        self.multiplesolutions = multiplesolutions
        self.allcyclesolutions = []
        self.solution_threshold = 150
        self.total_allowable_cyclecheck = 350
        self.variables_strings = [str(variable) for variable in LP.variables]
        optimalsolutions = []
        optimalsolutions_internal = []
   
        '''Set problem bounds and solve'''
        lp, obj = self.set_lp_problem(LP.lp, LP.variables)
        solution, solution_internal = self.ip_calculate(lp, LP.variables, obj)
        self.fill_allsolutions(solution)
        solution, solution_internal, lp, variables, obj =  self.initiate_cycle_check(solution, solution_internal, lp, LP.variables, obj, 0, initialcheck_value=True)
        optimalsolutions, optimalsolutions_internal, check_new_solution, check_solution_threshold =  self.filling_optimal_solution_arrays(solution, solution_internal, optimalsolutions, optimalsolutions_internal)
        optimalsolutions =  self.initiate_multiple_solutions(solution, lp, variables, obj, optimalsolutions, optimalsolutions_internal, 0)

        '''This is check that is needed only if multiple solutions were not identified and one wants to get then internal rxns of the one solution'''
        if multiplesolutions == 'False' and solution_internal:
            optimalsolutions, variables = self.identify_internal_rxns(LP.variables, optimalsolutions,
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
            verbose_print(self.verbose, 'Finding {} suboptimal shortest path for target {}'.format(count_k_paths, self.target))
            '''Solve for new pathway with new calculated weights'''

            for count_solution, solution in enumerate(self.allsolutions):
                if count_solution not in self.k_bounds:
                    temp = []
                    for r in solution:
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
            solution, solution_internal, lp, variables, obj =  self.initiate_cycle_check(solution, solution_internal, lp, variables, obj, 0, initialcheck_value=False)
            op, op_internal, check_new_solution, check_solution_threshold =  self.filling_optimal_solution_arrays(solution, solution_internal, op, op_internal)
            if solution and check_new_solution and check_solution_threshold:
                op = self.multiple_optimal_solution(lp, variables, obj, solution,
                                                    op, op_internal, count_k_paths)
            elif solution and (not check_new_solution or not check_solution_threshold):
                if op_internal:
                    op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
        return op
 
    def identify_internal_rxns(self, variables, op, op_internal, lp):
        internal_ignorerxns = []
        for orig_solution in op_internal:
            for variable in variables:
                rxn = self.allrxnsrev_dict[str(variable)]
                reaction = deepcopy(self.allrxnsrev_dict[str(variable)])
                reaction = re.sub('_F$', '', reaction)
                reaction = re.sub('_R$', '', reaction)
                if rxn in orig_solution:
                    variable.lowBound = 1
                    variable.upBound = 1
                elif reaction not in self.inrxns:
                    if variable.lowBound == 0 and variable.upBound == 0:
                        internal_ignorerxns.append(rxn) 
                    variable.lowBound = 0
                    variable.upBound = 0

            obj = self.set_objective_function_internal(variables)
            solution_orig, solution_internal = self.ip_calculate(lp, variables, obj)
            solution, lp, variables, obj = self.initiate_internal_cycle_check(solution_orig+solution_internal, lp, variables, obj, 0, len(solution_orig))
            op, op_internal, check_new_solution, check_solution_threshold = self.filling_optimal_solution_arrays(solution, [], op, op_internal)
            '''Add pathway to all optimal solutions'''

            if solution and check_solution_threshold and check_new_solution:
                verbose_print(self.verbose, 'STATUS: Checking for multiple optimal solutions internal ... ')
                op = self.multiple_optimal_solution_internal(lp, variables, obj, solution_orig, op)
            elif not check_solution_threshold:
                verbose_print(self.verbose, 'STATUS: Either solutions types internal {} and/or external {} have exceeded the solution threshold {} for target {}'.format(len(op_internal), len(op), self.solution_threshold, self.target))

            for variable in variables:
                rxn = self.allrxnsrev_dict[str(variable)]
                if rxn.startswith('EX_') or rxn in internal_ignorerxns:
                    variable.lowBound = 0
                    variable.upBound = 0

                else:
                    variable.lowBound = 0
                    variable.upBound = 1
        if not op:
            verbose_print(self.verbose, 'STATUS: No initial solution found when checking pathway with internal reactions for compound {} therefore trying again'.format(self.target))
            lp, obj = self.set_lp_problem(lp, variables)
            solution, solution_internal = self.ip_calculate(lp, variables, obj)
            self.fill_allsolutions(solution)
            solution, solution_internal, lp, variables, obj =  self.initiate_cycle_check(solution, solution_internal, lp, variables, obj, 0, initialcheck_value=True)
            op, op_internal, check_new_solution, check_solution_threshold =  self.filling_optimal_solution_arrays(solution, solution_internal, op, op_internal)
            op =  self.initiate_multiple_solutions(solution, lp, variables, obj, op, op_internal, 0)            
        return op, variables

    def multiple_optimal_solution_internal(self, lp, variables, obj, originalsolution, op):
        '''Find multiple solutions that use organisms internal reactions'''
        E = self.set_weight(len(originalsolution))
        for rxn in originalsolution:
            index = self.variables_strings.index(self.allrxnsrev_dict_rev[rxn])
            obj[index] = (1+E)
        solution_orig, solution_internal = self.ip_calculate(lp, variables, obj)
        if solution_orig not in op:
            solution, lp, variables, obj = self.initiate_internal_cycle_check(solution_orig+solution_internal, lp, variables, obj, 0, len(solution_orig))
            op, op_internal, check_new_solution, check_solution_threshold = self.filling_optimal_solution_arrays(solution, [], op, [])

            if len(solution_orig) == len(originalsolution):
                if check_new_solution and check_solution_threshold:
                    op = self.multiple_optimal_solution_internal(lp, variables, obj, solution, op)
                if not check_new_solution:
                    verbose_print(self.verbose, 'STATUS: solution with internal rxns for target {} not new therefore not adding to overall solutions ')
                if not check_solution_threshold:
                    verbose_print(self.verbose, 'STATUS: Either solutions types internal {} and/or external {} have exceeded the solution threshold {} for target {}'.format(len(op_internal), len(op), self.solution_threshold, self.target))
            else:
                verbose_print(self.verbose, 'STATUS: solution with internal rxns for target {} not the same length as original therefore not adding to overall solutions')
        else:
            verbose_print(self.verbose, 'STATUS: solution with internal rxns for target {} not new therefore not adding to overall solutions ')
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
        if solution not in op+op_internal:
            solution, solution_internal, lp, variables, obj = self.initiate_cycle_check(solution, solution_internal, lp, variables, obj, 0, initialcheck_value=False)
            op, op_internal, check_new_solution, check_solution_threshold = self.filling_optimal_solution_arrays(solution, solution_internal, op, op_internal)

            if len(solution) == len(originalsolution) and check_solution_threshold:
                '''If pathway not already identified add to total solution and check for more '''
                if check_new_solution:
                    op = self.multiple_optimal_solution(lp, variables, obj, solution, op, op_internal, count_k_paths)
                elif not check_new_solution:
                    if op_internal:
                        op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                    count_k_paths += 1
                    op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
            elif len(solution) != len(originalsolution) and check_solution_threshold:
                '''If solution is 0 or not the same length as previous solution check for next shortest solution'''
                if check_solution_threshold:
                    verbose_print(self.verbose, 'STATUS: No longer seeking multiple solutions for target {} at k{} level because solution {} is different length than original solution {}'.format(self.target, count_k_paths, solution, originalsolution))
                    if op_internal:
                        op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                    count_k_paths += 1
                    op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)

            elif len(solution) == len(originalsolution) and not check_solution_threshold:
                verbose_print(self.verbose, 'STATUS: Either solutions types internal {} and/or external {} have exceeded the solution threshold {} for target {}'.format(len(op_internal), len(op), self.solution_threshold, self.target))
                if op_internal:
                        op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)

            elif len(solution) != len(originalsolution) and not check_solution_threshold:
                verbose_print(self.verbose, 'STATUS: Either solutions types internal {} and/or external {} have exceeded the solution threshold {} for target {} or solutions do not match in length'.format(len(op_internal), len(op), self.solution_threshold, self.target))
                if op_internal:
                    op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
        else:
            if op_internal:
                op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
            count_k_paths += 1
            op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)          
        
        return op

    def cycle_constraints(self, lp, variables, solution, solution_internal, obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints and resolve
        if cycle is identified
        '''
        verbose_print(self.verbose, 'STATUS: Checking for cycles in the identified pathways for target {}'.format(self.target))

        '''Check if there is a cycle in identified pathway'''
        cycletest_1 = self.run_cycle_check(solution)
        original_solution = deepcopy(solution)

        if cycletest_1 is True:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if (self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle)) or cycle_count > self.total_allowable_cyclecheck:
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                if cycle_count > self.total_allowable_cyclecheck:
                    print ('STATUS: Exceeded number of code set cycle checks of {} for target {}'.format(self.total_allowable_cyclecheck, self.target))
                else:
                    print ('STATUS: Exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {} for target {}'.format(self.limit_cycle, self.target))
                return ([], [], lp, variables, obj)

            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''
                if solution not in self.allcyclesolutions:
                    temp = []
                    for r in solution:
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
                            print ('WARNING: new solution {} not the same number of reactions as the original {}  for target {}'.format(solution, original_solution, self.target))
                            return ([], [], lp, variables, obj)
                    else:
                        if solution:
                            verbose_print(self.verbose, 'STATUS: Check for cycle for solution {} for target {}'.format(solution, self.target))
                            solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables, solution,
                                                                                                     solution_internal, obj, cycle_count,
                                                                                                     initialcheck=True)
                            return (solution, solution_internal, lp, variables, obj)
                        else:
                            verbose_print(self.verbose, 'STATUS: No solution for target {}'.format(self.target))
                            return (solution, solution_internal, lp, variables, obj)

                else:
                    return ([], [], lp, variables, obj)

        elif cycletest_1 is False:
            return (solution, solution_internal, lp, variables, obj)

        elif cycletest_1 is None:
            print('WARNING: unable to perform cycle check for target {} consider decreasing processors'.format(self.target))
            return ([], [], lp, variables, obj)

    def cycle_constraints_internal(self, lp, variables, solution,  obj, cycle_count, length_external):
        '''
        Check solution for cycles and implement new constraints and resolve
        if cycle is identified
        '''
        verbose_print(self.verbose, 'STATUS: Checking for cycles in the identified internal pathways for target {}'.format(self.target))

        '''Check if there is a cycle in identified pathway'''
        cycletest_2 = self.run_cycle_check(solution)
        # original_solution = deepcopy(solution)

        if cycletest_2 is True:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if (self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle)) or cycle_count > self.total_allowable_cyclecheck:
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                if cycle_count > self.total_allowable_cyclecheck:
                    print ('STATUS: Exceeded number of code set cycle checks of {} for target {}'.format(self.total_allowable_cyclecheck, self.target))
                else:
                    print ('STATUS: Path (internal) exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {} for target {}'.format(self.limit_cycle, self.target))
                return ([], lp, variables, obj)

            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''
                if solution not in self.allcyclesolutions:
                    temp = []
                    for r in solution:
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append(self.variables_strings.index(self.allrxnsrev_dict_rev[r]))
                    self.allcyclesolutions.append(solution)
                    lp += pulp.LpConstraint(pulp.lpSum(1*variables[j] for j in temp), name='cycle pathway constraint '+str(len(self.allcyclesolutions)), sense=-1, rhs=len(temp)-1)
                    solution, solution_internal = self.ip_calculate(lp, variables, obj)
                    if len(solution) == length_external:
                        solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables,
                                                                                        solution+solution_internal,
                                                                                        obj, cycle_count, len(solution))
                        return (solution, lp, variables, obj)
                    else:
                        print ('WARNING: new solution {} different length than original {} for target {}'.format(len(solution), length_external, self.target))
                        return ([], lp, variables, obj)
        elif cycletest_2 is False:
            return(solution, lp, variables, obj)

        elif cycletest_2 is None:
            print('WARNING: unable to perform cycle check for target {} consider decreasing processors'.format(self.target))
            return ([], lp, variables, obj)

    def run_cycle_check(self, solution):
        '''
        Run cycle check
        '''
        if solution:
            totalvariables = []
            totalarcs = 0
            for rxn in solution:
                totalvariables.append(rxn)
                match = re.search('_F$', rxn)
                rxn_org = rxn
                if match is not None:
                    rxn_org = re.sub('_F$', '', rxn)
                match = re.search('_R$', rxn)
                if match is not None:
                    rxn_org = re.sub('_R$', '', rxn)
                reactants = self.DB.get_reactants(rxn_org)
                products = self.DB.get_products(rxn_org)
                if reactants is 'Errored' or products is 'Errored' or products is None or reactants is None:
                    print('WARNING: could not get reactants and/or products for cycle check for target {}'.format(self.target))
                    return (None)
                for reactant in reactants:
                    if reactant  not in self.incpds:
                        totalvariables.append(reactant)
                        totalarcs += 1
                for product in products:
                    if product not in self.incpds:
                        totalvariables.append(product)
                        totalarcs += 1
            # verbose_print(self.verbose, 'STATUS: for path {} arcs {}, variables {}, set variables {} for target {}'.format(solution, totalarcs, totalvariables, set(totalvariables), self.target))
            verbose_print(self.verbose, 'STATUS: optimal pathway {} for target {} has total arc count {} and total variable count -1 {}'.format(solution, self.target, totalarcs, len(set(totalvariables))-1))
            if totalarcs > len(set(totalvariables))-1:
                verbose_print(self.verbose, 'STATUS: optimal pathway {} has cycle for target {}'.format(solution, self.target))
                return (True)
            else:
                verbose_print(self.verbose, 'STATUS: No cycles were found in pathway {} for target {}'.format(solution, self.target))
                return (False)
        else:
            verbose_print(self.verbose, 'STATUS: pathway empty {} for target {}'.format(solution, self.target))
            return (False)