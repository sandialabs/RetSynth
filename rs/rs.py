from __future__ import print_function
__author__ = 'Leanne Whitmore and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov and cmhudso@sandia.gov'
__description__ = 'Main code to RetroSynthesis (RS)'

from multiprocessing import Process, Queue
import argparse
import cPickle
import os
import re
import glob
from Parser import read_startcompounds as rtsc
from Parser import read_targets as rt
from Parser import generate_output as go
from Visualization import SP_Graph_dot as spgd
from Visualization import reaction_files as rf
from ShortestPath import extractinfo as ei
from ShortestPath import constraints as co
from ShortestPath import integerprogram_glpk as ip_glpk
from ShortestPath import integerprogram_pulp as ip_pulp
from ShortestPath import search_sp_metclusters as smc
from Database import generate_database as gen_db
from Database import translate_metacyc as tm
from Database import query as Q
from FBA import build_model as bm
from FBA import optimize_target as ot
from FBA import compare_results as cr
from FBA import retrieve_producable_mets as rpm
from FBA import compareKO_results as crko
from RDFConverter import RDFileReaderMP
PATH = os.path.dirname(os.path.abspath(__file__))

def parse_arguments():
    '''
    The fundamental design of retrosynth is to take an organism, and a target
    chemical and output the minimum number of steps/reactions required to
    produce the target chemical in the organism.
    '''
    parser = argparse.ArgumentParser(description="BioRetroSynthesis Software: \
                                    Software identifies reactions and \
                                    corresponding genes that need to be added\
                                    to a desired organism to produce a target\
                                    chemical compound")
    parser.add_argument('-t', '--targets', help='Input file containing target \
                                                 compounds and organisms \
                                                 (tab deliminated)',
                        required=True, type=str)
    parser.add_argument('-gdb', '--generate_database', help='Generate database \
                        to use', type=str)
    parser.add_argument('-db', '--database', help='Specify database to use',
                        required=False, type=str)
    parser.add_argument('-mc', '--metacyc_addition', help='Add metacyc xml \
                        file to database', required=False, type=str)
    parser.add_argument('-rdf', '--rdf_addition', help='Add rdf folder to \
                        database', required=False, type=str)
    parser.add_argument('-tf', '--translation_file', help='Translation file to connect \
                        metacyc and kbase IDs', required=False, type=str)
    parser.add_argument('-d_dir', '--dump_directory', help='Path to folder \
                        of SBML network files(needed to generate database)',
                        required=False, type=str)
    parser.add_argument('-gdbc', '--generate_database_constraints',
                        help='Generate output constraint file for entire \
                        database', required=False, type=str)
    parser.add_argument('-dbc', '--database_constraints', help='Utilize \
                        constraint file for entire database',
                        required=False, type=str)
    parser.add_argument('-fba', '--flux_balance_analysis', help='Runs flux \
                        balance analysis using cobrapy on model organism of \
                        interest', required=False, action="store_true")
    parser.add_argument('-media', '--media_for_FBA', help='Define media for \
                        Flux Balance Analysis', required=False, type=str)
    parser.add_argument('-ko', '--knockouts', help='Performs single reaction \
                        knockouts on FBA model with added reaction and \
                        outputs rxns that result in an increase in the \
                        Objective Function', required=False,
                        action="store_true")
    parser.add_argument('--figures', help='Generate figures for clustering and \
                        filtering results', action='store_true')
    parser.add_argument('-lr', '--limit_reactions', help='Limit the number of \
                        reactions in a identified pathway (default:10, if no \
                        limit is wanted provide option of None)',
                        required=False, type=str, default=10)
    parser.add_argument('-lc', '--limit_cycles', help='Limits the number of cycle \
                        checks (default:10, if no limit is wanted provide \
                        option of None)  (not totally functional yet)', required=False, type=str,
                        default='None')
    parser.add_argument('--inchidb', help='Retrieve InChis and use them as compound \
                        IDs in the metabolic database', action='store_true')
    parser.add_argument('-op', '--output_path', help='Destination for output \
                        files', required=False, type=str, default=PATH)
    parser.add_argument('-rxntype', '--reaction_type', help='Define type of reactions \
                                                             being added to database \
                                                            (options are bio (default) and chem)',
                        required=False, type=str, default='bio')
    parser.add_argument('-rdfrxntype', '--rdf_reaction_type', help='Define type of reactions \
                                                                    from rdf files being \
                                                                    added to database (options \
                                                                    are bio and chem (default))',
                        required=False, type=str, default='chem')
    parser.add_argument('-evalrxns', '--evaluate_reactions', help='Defines which type of reactions \
                                                                   (bio, chem, or all (default)) \
                                                                   to be evaluated in identifying \
                                                                   shortest paths',
                        required=False, type=str, default='all')
    parser.add_argument('-stcpds', '--start_compounds', help='List of compounds that the shortest \
                                                              path should start from \
                                                              (use in place of target organism)',
                        required=False, type=str)
    parser.add_argument('-k', '--k_number_of_paths', help='Specifies the number of shortest \
                                                           paths wanted by the user',
                        required=False, type=int, default=0)
    parser.add_argument('-ms', '--multiple_solutions', help='Find all shortest paths (True \
                                                             (defualt) or False)',
                        required=False, type=str, default='True')
    parser.add_argument('-cy', '--cycles', help='Elimate cycles when finding shortest paths \
                                                 (True (default) or False)',
                        required=False, type=str, default='True')
    parser.add_argument('-p', '--processors', help='Number of processors to use when \
                                                    solving for shortest path (default 4)',
                        required=False, type=int, default=4)
    parser.add_argument('-solver', '--python_glpk_connector_package', help='python package to use \
                                                                            to connect to glpk \
                                                                            solver software (can \
                                                                            use GLPK (default) \
                                                                            or PULP)',
                        required=False, type=str, default='GLPK')
    parser.add_argument('--images', help='Set whether to use chemical images \
                                          or round nodes in output figures\
                                        (default True)', required=False,
                        type=str, default='True')
    return parser.parse_args()


def check_arguments(args):
    '''Checks and makes sure all required arguments are provided'''
    parser = argparse.ArgumentParser()
    if args.generate_database and not args.dump_directory:
        parser.error('If generating a new database, must specify a \
                      directory, --dump_directory, of metabolic networks \
                      (SBML)')
    if not args.generate_database_constraints \
            and not args.database_constraints:
        parser.error('Requires the use previously generated .constraints \
                     file for a database --database_constraints or generate a \
                     new .constraints file --generate_database_constraints')
    if not args.generate_database and not args.database:
        parser.error('Requires the use previously generated database \
                     --database or generate a database --generate_database')

    if args.metacyc_addition and not args.translation_file:
        parser.error('Requires the use --translation_file')

    if args.translation_file and not args.metacyc_addition:
        parser.error('Requires the use --metacyc_addition')

    if args.knockouts and not args.flux_balance_analysis:
        parser.error('--knockouts option requires that \
                     --flux_balance_analysis option be also specified')

    if args.media_for_FBA and not args.flux_balance_analysis:
        parser.error('--media_for_FBA option requires that \
                     --flux_balance_analysis option be also specified')

    if args.start_compounds and args.flux_balance_analysis:
        parser.error('Flux balance cannot be performed on a set of starting compounds\
                     would need to use an organisms metabolism to simulate flux')
    if not args.multiple_solutions and args.k_number_of_paths:
        parser.error('Cannot find k_number_of_paths correctly \
                     without finding all multiple_solutions')
    if not args.targets:
        parser.error('Requires an input file of target compounds')


def read_in_and_generate_output_files(args, database):
    '''Read in target input file and generate output files'''
    DB = Q.Connector(database)
    R = rt.Readfile(args.targets, DB, args.inchidb)
    if not R.targets:
        raise ValueError('ERROR: No targets, try different compounds')
    OUTPUT = go.Output(DB, args.output_path, args.flux_balance_analysis, args.knockouts)
    DB.conn.close()
    return(R.targets, R.ignorerxns, OUTPUT)

def retrieve_database_info(args):
    '''
    Generates database or uses previously generated database.
    Can also add metacyc database to a Kbase metabolic database
    '''
    if args.generate_database:
        '''
        Generate a database
        '''
        gen_db.Createdb(args.generate_database, args.dump_directory,
                        args.inchidb, args.reaction_type)
        DB = Q.Connector(args.generate_database)
        compartmentID = DB.get_compartment('cytosol')
        compartmentID = compartmentID[0]
        if args.metacyc_addition:
            '''
            Translate from metacyc database to database
            '''
            T = tm.Translate(args.generate_database, DB,
                             args.metacyc_addition, args.translation_file,
                             args.inchidb, args.reaction_type)
            DB = T.DB
        if args.rdf_addition:
            '''
            Translate synthetic rdf file to database
            '''
            RDFileReaderMP.RDF_Reader(args.rdf_addition,
                                      args.generate_database,
                                      args.rdf_reaction_type,
                                      compartmentID, args.processors)
            DB = Q.Connector(args.generate_database)
        database = args.generate_database
    elif args.database:
        '''
        Use and existing database
        '''
        DB = Q.Connector(args.database)
        compartmentID = DB.get_compartment('cytosol')
        compartmentID = compartmentID[0]
        if args.metacyc_addition:
            tm.Translate(args.database, DB, args.metacyc_addition,
                         args.translation_file, args.inchidb,
                         args.reaction_type)
        if args.rdf_addition:
            '''
            Add a synthetic rdf file to the database
            '''
            RDFileReaderMP.RDF_Reader(args.rdf_addition,
                                      args.database,
                                      args.rdf_reaction_type,
                                      compartmentID, args.processors)
            DB = Q.Connector(args.database)
        database = args.database

    allcpds = DB.get_all_compounds()
    if args.evaluate_reactions == 'all':
        allrxns = DB.get_all_reactions()
    elif args.evaluate_reactions == 'bio':
        allrxns = DB.get_reactions_based_on_type('bio')
    elif args.evaluate_reactions == 'chem':
        allrxns = DB.get_reactions_based_on_type('chem')
    DB.conn.close()
    return(allcpds, allrxns, database)

def retrieve_constraints(args, allrxns, allcpds, ignore_reactions, database):
    '''
    Generates database constraints or uses previously generated
    database constraints (.constraints) file
    '''
    DB = Q.Connector(database)
    if args.generate_database_constraints:
        LP = co.ConstructInitialLP(allrxns, allcpds, DB,
                                   ignore_reactions, True,
                                   reverse_constraints=False,
                                   specified_pysolver=args.python_glpk_connector_package)
        with open(args.generate_database_constraints, 'wb') as fout1:
            cPickle.dump(LP.A, fout1)
            cPickle.dump(LP.allcpds, fout1)

    elif args.database_constraints:
        with open(args.database_constraints, 'rb') as fin1:
            A = cPickle.load(fin1)
            allcompounds4matrix = cPickle.load(fin1)
        LP = co.ConstructInitialLP(allrxns, allcompounds4matrix, DB,
                                   ignore_reactions, A,
                                   reverse_constraints=False,
                                   specified_pysolver=args.python_glpk_connector_package)
    DB.conn.close()
    return LP

def construct_and_run_integerprogram(args, targets, LP, output, database):
    '''
    Constructs ILP and solves it identifying shortest path to the target
    '''
    DB = Q.Connector(database)
    if LP.PYSOLVER is 'GLPK':
        IP = ip_glpk.IntergerProgram(DB, args.limit_reactions,
                                     args.limit_cycles, args.k_number_of_paths, args.cycles)
    elif LP.PYSOLVER is 'PULP':
        IP = ip_pulp.IntergerProgram(DB, args.limit_reactions,
                                     args.limit_cycles, args.k_number_of_paths, args.cycles)
    else:
        raise IOError('ERROR: NO PYTHON SOLVER PACKAGE COULD BE IDENTIFIED. \
                      INSTALL PYGLPK, PYMPROG, OR PULP')
    if args.flux_balance_analysis:
        active_metabolism = retrieve_active_FBA_metabolism(targets, DB,
                                                           args, output)
    else:
        active_metabolism = {}
    return (IP, active_metabolism)

def _specific_target(target_id):
    '''Determines if there was a specified organism'''
    if target_id in ['', 'NA', 'N/A']:
        return False
    else:
        return True

def retrieve_shortestpath(target_info, IP, LP, database, args, output, active_metabolism):
    '''Retrieve the shortest path for target organism'''
    DB = Q.Connector(database)
    print (target_info)
    if args.images == 'False':
        _images = False
    else:
        _images = True
    if not _specific_target(target_info[2]) and not args.start_compounds:
        print ('WARNING: No organism given therefore software will run all organisms in database {} compound'.format(target_info[0]))
        ''' takes about 5 hours '''
        smc.SearchMetabolicClusters(target_info[0], LP, IP, output, DB)
    else:
        if args.start_compounds:
            incpds = rtsc.readfile_startcompounds(args.start_compounds)
            inrxns = []
        else:
            incpds = DB.get_compounds_in_model(target_info[2])
            inrxns = DB.get_reactions_in_model(target_info[2])

        if args.flux_balance_analysis:
            incpds_active = active_metabolism[target_info[2]][0]
            inrxns_active = active_metabolism[target_info[2]][1]
        else:
            incpds_active = incpds
            inrxns_active = inrxns

        if target_info[0] in incpds_active: #Check if compound exists in organism
            output.output_compound_natively_present_in_target_organism(target_info)
        else:
            optimal_pathways = IP.run_glpk(LP, incpds_active, inrxns_active, target_info[0],
                                           multiplesolutions=args.multiple_solutions)
            if optimal_pathways[0]:
                uniq_externalrxns = []
                for path in optimal_pathways:
                    path_org = []
                    for rxn in path:
                        rxn = re.sub('_F$', '', rxn)
                        rxn = re.sub('_R$', '', rxn)
                        path_org.append(rxn)
                    uniq_externalrxns.append(list(set(path_org) - set(inrxns)))

                ex_info = ei.Extract_Information(optimal_pathways, incpds_active, DB)
                output.output_shortest_paths(target_info, ex_info.temp_rxns)
                if args.flux_balance_analysis:
                    opt_fba = run_flux_balance_analysis(target_info, ex_info, incpds,
                                                        incpds_active, inrxns_active,
                                                        args.media_for_FBA, args.knockouts,
                                                        output, DB)
                    if args.figures:
                        G = spgd.GraphDot(DB, args.output_path, incpds_active, inrxns,
                                          opt_fba.fbasol.x_dict)
                        G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)
                        rf.ReactionFiles(args.output_path, DB, ex_info.temp_rxns,
                                         target_info[0], target_info[2], incpds_active)
                elif args.figures and not args.flux_balance_analysis:
                    G = spgd.GraphDot(DB, args.output_path, incpds, inrxns)
                    G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)
                    rf.ReactionFiles(args.output_path, DB, ex_info.temp_rxns,
                                     target_info[0], target_info[2], incpds)
            else:
                output.output_shortest_paths(target_info, optimal_pathways[0])
                if args.flux_balance_analysis:
                    print('WARNING: No optimal path for %s in species %s therefore no \
                          flux balance will be performed' % (target_info[0], target_info[2]))
    DB.conn.close()

def run_flux_balance_analysis(target_info, ex_info, incpds_original,
                              incpds_active, inrxns, media, ko,
                              output, DB):
    '''
    Run flux balance analysis on target organism with added reactions
    necessary to produce target compound
    '''
    fba = bm.BuildModel(target_info[2], incpds_original, inrxns, DB, media)
    opt_fba = ot.OptimizeTarget(target_info[0], target_info[2], fba.model, ex_info.temp_rxns,
                                ex_info.temp_exmets, fba.compounds_dict, incpds_active,
                                inrxns, DB, ko)
    comparisonresults = cr.Compare(target_info[0], fba.solution, opt_fba.fbasol,
                                   ex_info.temp_rxns, DB)
    output.output_FBA(target_info, fba.solution, opt_fba, comparisonresults, ex_info.temp_rxns)
    output.output_theoretical_yield(target_info[0], target_info[2], opt_fba.fbasol,
                                    opt_fba.compounds_dict)
    if ko:
        output.output_essential_reactions(target_info[0], target_info[2], opt_fba.essentialrxns)
        comparisonKOresults = crko.CompareKO(target_info[0], opt_fba.fbasol, opt_fba.KOsolutions,
                                             ex_info.temp_rxns, DB)
        output.output_FBA_KOs(target_info, comparisonKOresults, ex_info.temp_rxns)
    return opt_fba

def retrieve_active_FBA_metabolism(targets, DB, args, output):
    '''Retrieve active metabolism'''
    active_metabolism = {}
    organisms = set()
    output_queue = Queue()
    for target_info in targets:
        organisms.add(target_info[2])
    organisms = list(organisms)
    args_organisms = [organisms[i:i+args.processors]
                      for i in range(0, len(organisms), args.processors)]
    for orgs in args_organisms:
        processes = []
        for org in orgs:
            incpds = DB.get_compounds_in_model(org)
            inrxns = DB.get_reactions_in_model(org)
            processes.append(Process(target=rpm.RetrieveActiveRxnsCompounds,
                                     args=(org, incpds, inrxns, DB, output_queue,
                                           args.media_for_FBA)))
        for p in processes:
            p.start()
        for p in processes:
            ac = output_queue.get()
            active_metabolism.update(ac)
    output.output_activemetabolism(active_metabolism)
    return active_metabolism

def main():
    '''Main class'''
    args = parse_arguments()
    check_arguments(args)
    all_db_compounds, all_db_reactions, database = retrieve_database_info(args)
    targets, ignore_reactions, output = read_in_and_generate_output_files(args, database)
    LP = retrieve_constraints(args, all_db_reactions, all_db_compounds, ignore_reactions, database)
    IP, active_metabolism = construct_and_run_integerprogram(args, targets, LP, output, database)
    args_rcp_chunks = [targets[i:i+args.processors]
                       for i in range(0, len(targets), args.processors)]
    for chunks in args_rcp_chunks:
        processes = [Process(target=retrieve_shortestpath,
                             args=(chunk, IP, LP, database,
                                   args, output, active_metabolism)) for chunk in chunks]
        for p in processes:
            p.start()
        for p in processes:
            p.join()

    '''Remove all temporary images'''
    if args.images == 'True':
        for filename in glob.glob(PATH+"/Visualization/compound*"):
            os.remove(filename)

if __name__ == '__main__':
    main()
    exit()
