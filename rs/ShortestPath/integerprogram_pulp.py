from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Sets bounds necessary for a specific taraget compound for pulp and runs glpk'

import re
from copy import deepcopy
import pulp
from ShortestPath import cycle as cy
from ShortestPath import add_cycle_constraints_pulp as acc

class IntergerProgram(object):
    """Sets final constraints and solved integer linear program"""
    def __init__(self, db, limit_reactions, limit_cycle, k_paths, cycle, time_limit):
        '''initialize class'''
        self.limit_cycle = limit_cycle
        self.limit_reactions = limit_reactions
        self.k_paths = k_paths
        self.cycle = cycle
        self.time_limit = time_limit
        self.CYCLE = cy.CycleCheck(db)
        self.ACC = acc.AddCycleConstraints()
        self.weight_dict = {}

    def set_row_bounds(self, lp):
        '''Set row bounds'''
        print ('STATUS: Setting compound constraints ...')
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
        print ('STATUS: Generating objective function coefficients ...')
        obj = [0]*len(variables)
        for count, r in enumerate(self.allrxnsrev):
            reaction = deepcopy(r)
            reaction = re.sub('_F$', '', reaction)
            reaction = re.sub('_R$', '', reaction)
            if reaction not in self.inrxns:
                obj[count] = 1

        return obj
    def set_objective_function_internal(self, variables):
        '''Set objective function '''
        print ('STATUS: Generating objective function coefficients (internal)...')
        obj = [0]*len(variables)
        for count, r in enumerate(self.allrxnsrev):
            obj[count] = 1
        return obj

    def fill_allsolutions(self, solution):
        if solution not in self.allsolutions:
            self.allsolutions.append(solution)

    def ip_calculate(self, lp, variables, obj):
        '''Run glpk solver'''
        solution = []
        solution_internal = []
        print ('STATUS: Setting objective function ...')
        lp.setObjective(pulp.lpSum(obj[i]*variables[i] for i in range(len(obj))))

        print ('STATUS: Solving problem ...')
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
                    rxn_sol = self.allrxnsrev[self.LP.rxnnames.index(variable.name)]
                    if rxn_sol.endswith('_F'):
                        rxn_sol = re.sub('_F$', '', rxn_sol)
                    elif rxn_sol.endswith('_R'):
                        rxn_sol = re.sub('_R$', '', rxn_sol)
                    if rxn_sol not in self.inrxns:
                        solution.append(self.allrxnsrev[self.LP.rxnnames.index(variable.name)])
                    else:
                        solution_internal.append(self.allrxnsrev[self.LP.rxnnames.index(variable.name)])
        if self.limit_reactions != 'None' and len(solution) > int(self.limit_reactions):
            print ('STATUS: Path contains too many reaction steps, consider increasing limit up from {}'.format(self.limit_reactions))
            solution = []
            solution_internal = []
        return solution, solution_internal

    def run_glpk(self, LP, incpds, inrxns, target_compound_ID, multiplesolutions=True):
        '''Final set up and solve integer linear program'''
        '''Set initial variables'''
        self.inrxns = inrxns
        self.incpds = incpds
        self.target = target_compound_ID
        self.LP = LP
        self.allrxnsrev = self.LP.allrxnsrev
        self.allcpds = self.LP.allcpds
        lp = self.LP.lp
        self.allsolutions = []
        self.k_bounds = []
        optimalsolutions = []
        optimalsolutions_internal = []


        '''Set problem bounds and solve'''
        lp = self.set_row_bounds(lp)
        obj = self.set_objective_function(self.LP.variables)
        solution, solution_internal = self.ip_calculate(lp, self.LP.variables, obj)
        self.fill_allsolutions(solution)
        if self.cycle == 'True' and solution:
            '''Check for cycles in pathway'''
            solution, solution_internal, lp, self.LP.variables, obj = self.cycle_constraints(lp, self.LP.variables,
                                                                          solution, solution_internal, obj, 0,
                                                                          initialcheck=True)
        if solution_internal:
            self.fill_allsolutions(solution)
            optimalsolutions_internal.append(solution)
        else:
            self.fill_allsolutions(solution)
            optimalsolutions.append(solution)            

        if multiplesolutions == 'True' and solution:
            '''Check for multiple solutions'''
            print ('STATUS: Checking for multiple optimal solutions ... ')
            optimalsolutions = self.multiple_optimal_solution(lp, self.LP.variables,
                                                              obj, solution, optimalsolutions, optimalsolutions_internal, 0)
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
            print ('Finding {} shortest path'.format(count_k_paths))
            '''Solve for new pathway with new calculated weights'''
            
            for count_solution, solution in enumerate(self.allsolutions):
                if count_solution not in self.k_bounds:
                    temp = []
                    for count, r in enumerate(solution):
                        reaction = deepcopy(r)
                        reaction = re.sub('_F$', '', reaction)
                        reaction = re.sub('_R$', '', reaction)
                        if reaction not in self.inrxns:
                            temp.append(self.allrxnsrev.index(r))
                    lp += pulp.LpConstraint(pulp.lpSum(1*variables[j] for j in temp), name='K pathway constraint '+str(count_solution), sense=-1, rhs=len(temp)-1)
                    self.k_bounds.append(count_solution)
            solution, solution_internal = self.ip_calculate(lp, variables, obj)
            self.fill_allsolutions(solution)

            '''If pathway is greater in steps than the previous path add too optimal solutions
                check for multiple solution'''
            if self.cycle == 'True':
                solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables,
                                                                                         solution, solution_internal, obj, 0)
            if solution_internal:
                if solution and solution not in op_internal:
                    op_internal.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj, solution, op, op_internal, count_k_paths)
                else:
                    op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
            else:
                if solution:
                    op.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj,
                                                        solution, op, op_internal, count_k_paths)
        return op
 
    def identify_internal_rxns(self, variables, op, op_internal, lp):
        print ('IDENTIFY INTERNAL REACTIONS')
        ignorerxns = []
        for orig_solution in op_internal:
            print (str(orig_solution)+' op internal solution')
            for count, rxn in enumerate(self.allrxnsrev):
                reaction = deepcopy(rxn)
                reaction = re.sub('_F$', '', reaction)
                reaction = re.sub('_R$', '', reaction)
                if rxn in orig_solution:
                    variables[count].lowBound = 1
                    variables[count].upBound = 1
                elif reaction not in self.inrxns:
                    if variables[count].lowBound == 0 and variables[count].upBound == 0:
                        ignorerxns.append(rxn) 
                    variables[count].lowBound = 0
                    variables[count].upBound = 0

            obj = self.set_objective_function_internal(variables)
            solution, solution_internal = self.ip_calculate(lp, variables, obj)
            if self.cycle == 'True' and solution:
                '''Check for cycles in pathway'''
                solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables, solution+solution_internal,
                                                                                                  obj, 0, initialcheck=True)
                op.append(solution)
            else:
                solution = solution + solution_internal
                op.append(solution)

            '''Add pathway to all optimal solutions'''

            if solution:
                print ('STATUS: Checking for multiple optimal solutions internal ... ')
                op = self.multiple_optimal_solution_internal(lp, variables, obj, solution, op)
                print (op)

            for count, rxn in enumerate(self.allrxnsrev):
                if rxn.startswith('EX_') or rxn in ignorerxns:
                    variables[count].lowBound = 0
                    variables[count].upBound = 0

                else:
                    variables[count].lowBound = 0
                    variables[count].upBound = 1
        return op, variables

    def multiple_optimal_solution_internal(self, lp, variables, obj, originalsolution, op):
        '''Find multiple solutions that use organisms internal reactions'''
        E = self.set_weight(len(originalsolution))
        for rxn in originalsolution:
            index = self.allrxnsrev.index(rxn)
            obj[index] = (1+E)
        solution, solution_internal = self.ip_calculate(lp, variables, obj)
        if self.cycle =='True' and solution_internal:
            solution, lp, variables, obj = self.cycle_constraints_internal(lp, variables, solution+solution_internal, obj, 0)
        else:
            solution = solution+solution_internal
        if len(solution) == len(originalsolution):
            if solution and solution not in op:
                op.append(solution)
                op = self.multiple_optimal_solution_internal(lp, variables, obj, solution, op)
        return(op)

    def multiple_optimal_solution(self, lp, variables, obj, originalsolution, op, op_internal, count_k_paths):
        '''Identify multiple solutions'''
        '''Get weight'''
        E = self.set_weight(len(originalsolution))

        '''Set weight for reactions in objective function'''
        for rxn in originalsolution:
            index = self.allrxnsrev.index(rxn)
            obj[index] = (1+E)
        solution, solution_internal = self.ip_calculate(lp, variables, obj)
        self.fill_allsolutions(solution)
        '''Check for cycles'''
        if self.cycle == 'True':
            solution, solution_internal, lp, variables, obj = self.cycle_constraints(lp, variables,
                                                                                    solution, solution_internal, obj, 0)

        if len(solution) == len(originalsolution):
            '''If pathway not already identified add to total solution and check for more '''
            if solution_internal:
                if solution and solution not in op_internal:
                    op_internal.append(solution)
                    self.fill_allsolutions(solution)
                    op = self.multiple_optimal_solution(lp, variables, obj, solution, op, op_internal, count_k_paths)
                else:
                   op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                   count_k_paths += 1
                   op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
            else:
                if solution and solution+solution_internal not in op:
                    op.append(solution+solution_internal)
                    self.fill_allsolutions(solution)
                    if op_internal:
                        op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                    op = self.multiple_optimal_solution(lp, variables, obj,
                                                        solution, op, [], count_k_paths)

                else:
                    '''If pathway already in solution check for next shortest pathway'''
                    if op_internal:
                        op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
                    count_k_paths += 1
                    op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
        else:
            '''If solution is 0 or not the same length as previous solution check for next
                shortest solution'''
            if op_internal:
                op, variables = self.identify_internal_rxns(variables, op, op_internal, lp)
            count_k_paths += 1
            op = self.k_number_paths(lp, variables, obj, op, [], count_k_paths)
            print ('STATUS: No more multiple solutions...')
        return op

    def cycle_constraints(self, lp, variables, solution, solution_internal, obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints and resolve
        if cycle is identified
        '''
        print ('STATUS: Checking for cycles in the identified pathways')

        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ('STATUS: Exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {}'.format(self.limit_cycle))
                return ([], [], lp, variables, obj)

            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''

                print ('STATUS: Implementing new constraints for pathways with cycles')
                self.ACC.add_cycle_constraints(lp, variables, self.allrxnsrev, self.CYCLE)

                print ('STATUS: Solving new integer program with cycle constraints')
                add2obj = len(self.ACC.variables)-len(obj)
                add2obj_list = [0]*add2obj
                obj = obj + add2obj_list
                solution, solution_internal = self.ip_calculate(self.ACC.lp, self.ACC.variables, obj)

                '''Checks and makes sure resulting solution is of equal length as previous'''
                if initialcheck is False:
                    if len(solution) == len(original_solution):
                        solution, solution_internal, self.ACC.lp, self.ACC.variables, obj = self.cycle_constraints(self.ACC.lp, self.ACC.variables, solution, solution_internal, obj, cycle_count)
                        return (solution, solution_internal, self.ACC.lp, self.ACC.variables, obj)
                    else:
                        print ('WARNING: new solution longer than original {}'.format(solution))
                        return ([], [], self.ACC.lp, self.ACC.variables, obj)
                else:
                    solution, solution_internal, self.ACC.lp, self.ACC.variables, obj = self.cycle_constraints(self.ACC.lp, self.ACC.variables, solution, solution_internal, obj, cycle_count, initialcheck=True)
                    return (solution, solution_internal, self.ACC.lp, self.ACC.variables, obj)
        else:
            return (solution, solution_internal, lp, variables, obj)

    def cycle_constraints_internal(self, lp, variables, solution,  obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints and resolve
        if cycle is identified
        '''
        print ('STATUS: Checking for cycles in the identified pathways')

        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ('STATUS: Exceeded number of cycle checks, no pathways without cycles, consider increasing limit up from {}'.format(self.limit_cycle))
                return ([], lp, variables, obj)

            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''

                print ('STATUS: Implementing new constraints for pathways with cycles')
                self.ACC.add_cycle_constraints(lp, variables, self.allrxnsrev, self.CYCLE)

                print ('STATUS: Solving new integer program with cycle constraints')
                add2obj = len(self.ACC.variables)-len(obj)
                add2obj_list = [0]*add2obj
                obj = obj + add2obj_list
                solution, solution_internal = self.ip_calculate(self.ACC.lp, self.ACC.variables, obj)

                '''Checks and makes sure resulting solution is of equal length as previous'''
                if initialcheck is False:
                    if len(solution+solution_internal) == len(original_solution):
                        solution, self.ACC.lp, self.ACC.variables, obj = self.cycle_constraints_internal(self.ACC.lp, self.ACC.variables, solution+solution_internal, obj, cycle_count)
                        return (solution, self.ACC.lp, self.ACC.variables, obj)
                    else:
                        print ('WARNING: new solution longer than original {}'.format(solution))
                        return ([], self.ACC.lp, self.ACC.variables, obj)
                else:
                    solution, self.ACC.lp, self.ACC.variables, obj = self.cycle_constraints_internal(self.ACC.lp, self.ACC.variables, solution+solution_internal, obj, cycle_count, initialcheck=True)
                    return (solution, self.ACC.lp, self.ACC.variables, obj)
        else:
            return (solution, lp, variables, obj)

