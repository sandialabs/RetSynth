from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Sets bounds necessary for a specific taraget compound for pyglpk and runs glpk'

from copy import deepcopy
import re
from ShortestPath import cycle as cy
from ShortestPath import add_cycle_constraints_glpk as acc

class IntergerProgram(object):
    """Sets final constraints and solved integer linear program"""
    def __init__(self, db, limit_reactions, limit_cycle, k_paths, cycle):
        '''initalize class'''
        self.limit_cycle = limit_cycle
        self.limit_reactions = limit_reactions
        self.k_paths = k_paths
        self.cycle = cycle
        self.CYCLE = cy.CycleCheck(db)
        self.ACC = acc.AddCycleConstraintsGLPK()
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

    def ip_calculate(self, lp, obj):
        '''Run glpk solver'''
        solution = []
        lp.obj[:] = obj
        lp.intopt()
        if lp.status == 'undef':
            '''Check if solution was found'''
            pass
        else:
            '''If it was get reactions that were found to be in the shortest path'''
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
            if self.limit_reactions != 'None' and len(solution) > int(self.limit_reactions):
                print ('STATUS: Path contains too many reaction steps, consider increasing limit up from {} first line'.format(self.limit_reactions))
                solution = []
        return solution

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
        lp = self.LP.lp
        optimalsolutions = []

        '''Set problem bounds and solve'''
        lp = self.set_row_bounds(lp)
        obj = self.set_objective_function(lp)
        solution = self.ip_calculate(lp, obj)

        if self.cycle == 'True' and solution:
            '''Check for cycles in pathway'''
            solution, lp, obj = self.cycle_constraints(lp, solution, obj, 0, initialcheck=True)

        '''Add pathway to all optimal solutions'''
        optimalsolutions.append(solution)

        if multiplesolutions == 'True' and solution:
            '''Check for multiple solutions'''
            print ('STATUS: Checking for multiple optimal solutions ... ')
            optimalsolutions = self.multiple_optimal_solution(lp, obj, solution,
                                                              optimalsolutions, 0)

        return optimalsolutions

    def set_weight(self, number_rxn_steps):
        '''retrieve weight to disfavor rxns in path'''
        E = float(1)/float(number_rxn_steps)
        E = round(E, 2)
        E = E-.01
        return E

    def set_weights_k_paths(self, obj, new_weight, original_weight, original_num_steps):
        '''set weights to find the next shortest path '''
        '''Add new weight to weight dictionary'''
        self.weight_dict.setdefault(original_num_steps, original_weight+new_weight)

        for i, value in enumerate(obj):
            if value == original_weight:
                '''If value in obj equals the original weight add the new weight
                    too the value'''
                obj[i] = obj[i]+new_weight

            elif value in self.weight_dict.values():
                '''If a value in obj function equals a weight in the
                    weight dictionary add new weight plus original weight
                    to obj value'''
                obj[i] = obj[i]+original_weight+new_weight

        '''Add new weights too dictionary'''
        for key in self.weight_dict:
            if key != original_num_steps:
                self.weight_dict[key] = self.weight_dict[key]+(original_weight+new_weight)

        return obj

    def check_for_next_length_k_paths(self, obj, new_weight, original_weight,
                                      lp, count_k_paths, add_length, op):
        '''Executed if no pathway could be identified and the first given length'''
        if (len(op[-1])+add_length) > self.limit_reactions:
            '''Check to see if pathways identifed would exceed reaction limit,
                if so, stop looking for pathways'''
            print ('STATUS: Any pathways identified would be greater than reaction limit of {}, stop looking for solutions'.format(self.limit_reactions))
        else:
            print ('STATUS: No pathways of the length {} looking for pathways of length'.format(len(op[-1])+add_length),
                   len(op[-1])+1+add_length)
            E = self.set_weight(len(op[-1])+add_length)
            original_weight = 1+E

            '''set new weights for next shortest path'''
            obj = self.set_weights_k_paths(obj, new_weight, original_weight,
                                           len(op[-1])+add_length)

            '''Find next shortest path'''
            op = self.k_number_paths(lp, obj, op, original_weight,
                                     count_k_paths, add_length=add_length)
        return op

    def k_number_paths(self, lp, obj, op, original_weight, count_k_paths, add_length=0):
        ''' 'retrieve the next shortest path '''
        if count_k_paths <= self.k_paths:
            '''Calculate new weight'''
            new_weight = (float(1)/float(len(op[-1])+1)+add_length)+.02
            '''Apply new weights to objective function'''
            obj = self.set_weights_k_paths(obj, new_weight, original_weight, len(op[-1]))
            '''Solve for new pathway with new calculated weights'''
            solution = self.ip_calculate(lp, obj)
            if len(solution) <= len(op[-1]):
                '''If pathway is less or the same length as previously identified pathway get
                    weight for the next shortest path and resolve'''
                add_length += 1
                op = self.check_for_next_length_k_paths(obj, new_weight, original_weight,
                                                        lp, count_k_paths, add_length, op)

            else:
                '''If pathway is greater in steps than the previous path,
                 check for cycles and then add too optimal solutions
                 check for multiple solution'''
                if self.cycle == 'True':
                    solution, lp, obj = self.cycle_constraints(lp, solution, obj, 0)

                if solution:
                    op.append(solution)
                    op = self.multiple_optimal_solution(lp, obj, solution, op, count_k_paths)
                else:
                    print ('STATUS:  Cycles were identified in all pathways of this length, getting pathways of the next length')
                    add_length += 1
                    op = self.check_for_next_length_k_paths(add_length, op)

        return op

    def multiple_optimal_solution(self, lp, obj, originalsolution, op, count_k_paths):
        '''Identify multiple solutions'''
        '''Get weight'''
        E = self.set_weight(len(originalsolution))

        '''Set weight for reactions in objective function'''
        for rxn in originalsolution:
            index = self.allrxnsrev.index(rxn)
            obj[index] = (1+E)
        solution = self.ip_calculate(lp, obj)

        '''Check for cycles'''
        if self.cycle == 'True':
            solution, lp, obj = self.cycle_constraints(lp, solution, obj, 0)

        if len(solution) == len(originalsolution):
            ''' Check if solution is of the same number of reaction steps and does
                not have 0 reaction steps'''

            if solution and solution not in op:
                '''If pathway not already identified add to total solution and check for more '''
                op.append(solution)
                op = self.multiple_optimal_solution(lp, obj, solution, op, count_k_paths)

            else:
                '''If pathway already in solution check for next shortest pathway'''
                count_k_paths += 1
                op = self.k_number_paths(lp, obj, op, (1+E), count_k_paths)

        else:
            '''If solution is 0 or not the same length as previous solution check for next
                shortest solution'''
            count_k_paths += 1
            op = self.k_number_paths(lp, obj, op, (1+E), count_k_paths)
            print ('STATUS: No more multiple solutions...')
        return op

    def cycle_constraints(self, lp, solution, obj, cycle_count, initialcheck=False):
        '''
        Check solution for cycles and implement new constraints
        and resolve if cycle is identified
        '''
        print ('STATUS: Checking for cycles in the identified pathway')

        '''Check if there is a cycle in identified pathway'''
        cycletest = self.CYCLE.run_cycle_check(solution, self.incpds)
        original_solution = deepcopy(solution)

        if cycletest:
            '''If pathway has cycle begin cycle elimination steps'''
            cycle_count += 1
            if self.limit_cycle != 'None' and cycle_count > int(self.limit_cycle):
                '''Count number of cycle checks performed, if it exceeds designated limit
                    of cycle checks stop check for cycles and returns no path'''
                print ("STATUS: Exceeded number of cycle checks, \
                       no pathways without cycles, consider increasing limit up from {} ".format(self.limit_cycle))
                return ([], lp, obj)
            else:
                '''If cycle checks have not been exceed or have been set to None set new constraints
                    eliminate identification of pathway and solve for new shortest path'''
                print ('STATUS: Implementing new constraints for pathways with cycles')
                self.ACC.add_cycle_constraints(lp, self.allrxnsrev, self.CYCLE)

                print ('STATUS: Solving new integer program with cycle constraints')
                add2obj = len(self.ACC.lp.cols)-len(obj)
                add2obj_list = [0]*add2obj
                obj = obj + add2obj_list
                solution = self.ip_calculate(self.ACC.lp, obj)
                '''Checks and makes sure resulting solution is of equal length as previous'''
                if initialcheck is False:
                    if len(solution) == len(original_solution):
                        solution, self.ACC.lp, obj = self.cycle_constraints(self.ACC.lp,
                                                                            solution, obj,
                                                                            cycle_count)
                        return (solution, self.ACC.lp, obj)
                    else:
                        print ('WARNING: new solution longer than original {}'.format(solution))
                        return ([], lp, obj)
                else:
                    solution, self.ACC.lp, obj = self.cycle_constraints(self.ACC.lp, solution,
                                                                        obj, cycle_count,
                                                                        initialcheck=True)
                    return (solution, self.ACC.lp, obj)

        else:
            return (solution, lp, obj)
