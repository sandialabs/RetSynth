from __future__ import print_function
__author__ = 'Leanne Whitmore and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov and cmhudso@sandia.gov'
__description__ = 'Main code to RetSynth (RS)'

from multiprocessing import Process, Queue
import argparse
import cPickle
import os
import re
import glob
import time
from Parser import read_startcompounds as rtsc
from Parser import read_targets as rt
from Parser import generate_output as go
from Parser import structure_similarity as ss
from Visualization import SP_Graph_dot as spgd
from Visualization import reaction_files as rf
from ShortestPath import extractinfo as ei
from ShortestPath import constraints as co
from ShortestPath import integerprogram_pulp as ip_pulp
from ShortestPath import search_sp_metclusters as smc
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import build_KEGG_db as bkeggdb
from Database import build_metacyc_db as bmcdb
from Database import build_ATLAS_db as batlasdb
from Database import build_MINE_db as bminedb
from Database import build_SPRESI_db as bspresidb
from Database import query as Q
from Database import remove_duplicate_cpds as rdc
from FBA import build_model as bm
from FBA import optimize_target as ot
from FBA import compare_results as cr
from FBA import retrieve_producable_mets as rpm
from FBA import compareKO_results as crko
#from toxicity import TrainToxModel as tt

PATH = os.path.dirname(os.path.abspath(__file__))

def verbose_print(verbose, line):
    if verbose:
        print(line)

def parse_arguments():
    '''
    The fundamental design of retrosynth is to take an organism, and a target
    chemical and output the minimum number of steps/reactions required to
    produce the target chemical in the organism.
    '''
    parser = argparse.ArgumentParser(description="RetSynth Software: \
                                    Software identifies reactions and \
                                    corresponding genes that need to be added\
                                    to a desired organism to produce a target\
                                    chemical compound")
    parser.add_argument('-t', '--targets', help='Input file containing target \
                                                 compounds and organisms \
                                                 (tab deliminated)',
                        required=True, type=str)
    ###SOLVER OPTIONS### 
    parser.add_argument('-tmlim', '--solver_time_limit', help='time limit for solver to solve shortest \
                                                                path (note: only function with PULP python \
                                                                solver package)',
                                                                required=False, type=str, default=30)

    ###DATABASE OPTIONS###
    parser.add_argument('-gdb', '--generate_database', help='Generate database \
                        to use', type=str)
    parser.add_argument('-db', '--database', help='Specify database to use',
                        required=False, type=str)
    parser.add_argument('-dbc', '--database_constraints', help='Utilize \
                        constraint file for entire database',
                        required=False, type=str)    
    parser.add_argument('-gdbc', '--generate_database_constraints',
                        help='Generate output constraint file for entire \
                        database', required=False, type=str)

    parser.add_argument('--kbase', help='Set whether to build database with  Kbase data',  required=False,
                        action="store_true")
    
    parser.add_argument('-k_dir', '--kbase_dump_directory', help='Path to folder \
                        of SBML network files from kbase',
                        required=False, type=str)

    parser.add_argument('--metacyc', help='Set whether to build database with metacyc',  required=False,
                        action="store_true")    
    
    parser.add_argument('-mc', '--metacyc_addition', help='Add metacyc xml \
                        file to database', required=False, type=str)
    
    parser.add_argument('-tf', '--translation_file', help='Translation file to connect \
                        metacyc and kbase IDs', required=False, type=str)    


    parser.add_argument('--SPRESI', help='Set whether to build database with SPRESI information',
                        required=False, action="store_true")
    
    parser.add_argument('-s_dir', '--spresi_dump_directory', help='Path to folder of \
                        SPRESI files (.rdf)', required=False, type=str)

    parser.add_argument('--mine', help='Set whether to build database with mine data',  required=False,
                        action="store_true")
    
    parser.add_argument('-m_dir', '--mine_dump_directory', help='Path to folder \
                        of .msp files from MINE database',
                        required=False, type=str)

    parser.add_argument('--atlas', help='Set whether to build database with atlas data',  required=False,
                        action="store_true")
    parser.add_argument('-a_dir', '--atlas_dump_directory', help='Path to folder \
                        of ATLAS files (.csv)', required=False, type=str)

    parser.add_argument('--kegg', help='Set whether to build database with Kegg data',  required=False,
                        action="store_true")

    parser.add_argument('-kbaserxntype', '--kbase_reaction_type', help='Define type of reactions \
                                                             that are being added to database \
                                                            (options are bio (default) and chem)',
                        required=False, type=str, default='bio')
    parser.add_argument('-mcrxntype', '--metacyc_reaction_type', help='Define type of reactions \
                                                             that are being added to database \
                                                            (options are bio (default) and chem)',
                        required=False, type=str, default='bio')

    parser.add_argument('-keggrxntype', '--kegg_reaction_type', help='Define type of reactions \
                                                             that are being added to database \
                                                            (options are bio (default) and chem)',
                        required=False, type=str, default='bio')

    parser.add_argument('-keggorganismtype', '--kegg_organism_type', help='Define type of organisms \
                                                                           reactions to be added to \
                                                                           the database bacteria (default), \
                                                                           algae, plants, fungi, Eukaryotes, \
                                                                           prokaryotes, or all.  If user wants both plants \
                                                                           and bacteria seperate list by comma i.e. plants,bacteria', 
                        required=False, type=str, default='bacteria')

    parser.add_argument('-keggnunorganisms', '--kegg_number_of_organisms', help='Define number of organisms \
                                                                                    from kegg database to add',
                        required=False, type=str, default='all')

    parser.add_argument('-keggnunorganismpaths', '--kegg_number_of_organism_pathways', help='Define number of pathways \
                                                                                             from an organism to add',
                        required=False, type=str, default='all')

    parser.add_argument('-spresirxntype', '--spresi_reaction_type', help='Define type of reactions \
                                                                    from rdf files that are being \
                                                                    added to database (options \
                                                                    are bio and chem (default))',
                        required=False, type=str, default='chem')

    parser.add_argument('-atlasrxntype', '--atlas_reaction_type', help='Define type of reactions \
                                                                    from atlas files that are being \
                                                                    added to database (options \
                                                                    are bio (default) and chem)',
                        required=False, type=str, default='bio')

    parser.add_argument('-minerxntype', '--mine_reaction_type', help='Define type of reactions \
                                                                    from atlas files that are being \
                                                                    added to database (options \
                                                                    are bio (default) and chem)',
                        required=False, type=str, default='bio')

    ###FLUX BALANCE ANALYSIS OPTIONS###    
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

    ###EXTRA SOLVER OPTIONS###
    parser.add_argument('-lr', '--limit_reactions', help='Limit the number of reactions in a\
                                                          identified pathway (default:10, if no \
                                                          limit is wanted provide option of None)',
                        required=False, type=str, default=10)

    parser.add_argument('-lc', '--limit_cycles', help='Limits the number of cycle \
                                                       checks (default:10, if no limit is wanted provide \
                                                      option of None)  (not totally functional yet)',
                        required=False, type=str, default='None')

    parser.add_argument('--inchidb', help='Retrieve InChis and use them as compound \
                        IDs in the metabolic database', action='store_true')

    parser.add_argument('-op', '--output_path', help='Destination for output files',
                        required=False, type=str, default=PATH)

    parser.add_argument('-oxf', '--output_xlsx_format', help='Convert output txt files to xlsx \
                                                              files (warning: do not do this if \
                                                              getting pathways for 100s of compounds',
                        required=False, action="store_true")

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

    parser.add_argument('-tan_thresh', '--tanimoto_threshold', help='Set the tanimoto threshold to \
                                                                     identify related compounds',
                         required=False, type=str, default=1)

    parser.add_argument('-p', '--processors', help='Number of processors to use when \
                                                    solving for shortest path (default 4)',
                        required=False, type=int, default=4)

    ###FIGURE OPTIONS###
    parser.add_argument('--figures', help='Generate figures for clustering and \
                        filtering results', action='store_true')

    parser.add_argument('--images', help='Set whether to use chemical images \
                                          or round nodes in output figures\
                                          (default True)', required=False,
                        type=str, default='True')

    ###TOXICITY OPTIONS###
    parser.add_argument('-toxicity','--predict_toxicity',help='Predict whether non native compounds are \
                                                               toxic to host organism (currently only works with E. Coli Strains)',
                        required=False, action="store_true")
    parser.add_argument('-v', '--verbose', help="Level of output", action="store_true")    
    return parser.parse_args()


def check_arguments(args):
    '''Checks and makes sure all required arguments are provided'''
    parser = argparse.ArgumentParser()
    if not args.generate_database_constraints \
            and not args.database_constraints:
        parser.error('Requires the use previously generated .constraints \
                     file for a database --database_constraints or generate a \
                     new .constraints file --generate_database_constraints')
    if not args.generate_database and not args.database:
        parser.error('Requires the use previously generated database \
                     --database or generate a database --generate_database')

    if  args.generate_database and not (args.kbase or args.metacyc or args.kegg):
        parser.error('Requires a specified database type --kbase, --metacyc or --kegg')

    if args.kbase and not args.kbase_dump_directory:
        parser.error('Requires use of --kbase_dump_directory')

    if args.atlas and not args.atlas_dump_directory:
        parser.error('Requires use of --atlas_dump_directory')

    if args.SPRESI and not args.spresi_dump_directory:
        parser.error('Requires use of --spresi_dump_directory')

    if args.mine and not args.mine_dump_directory:
        parser.error('Requires use of --mine_dump_directory')

    if args.metacyc and not args.translation_file and not args.metacyc_addition:
        parser.error('Requires the use of --translation_file and --metacyc_addition')

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


def get_new_temp_imgs_folder(count):
    '''Check if folder to store images is already present if so new temp folder is made'''
    count+=1
    try:
        os.mkdir(PATH+'/temp_imgs_'+str(count))
        return PATH+'/temp_imgs_'+str(count)
    except OSError:
        PATH_NEW = get_new_temp_imgs_folder(count)
        return PATH_NEW
def get_compartmentID_from_db(DB, compartment):
    '''Retrieves specified compartment ID'''
    compartment = compartment.lower()
    compartmentID_array = DB.get_compartment(compartment)
    if compartmentID_array is None or len(compartmentID_array) == 0 or compartmentID_array[0] == '':
        print ('WARNING: Could not retrieve a compartment ID from the database')
        if compartment == 'cytosol':
            compartmentID = 'c0'
        elif compartment == 'extracellular':
            compartmentID = 'e0'
        else:
            compartmentID = 'c0'
    else:
        compartmentID = compartmentID_array[0]
    return (compartmentID)

def read_in_and_generate_output_files(args, database):
    '''Read in target input file and generate output files'''
    DB = Q.Connector(database)
    R = rt.Readfile(args.targets, DB, args.inchidb)
    if not R.targets:
        raise ValueError('ERROR: No targets, try different compounds')
    temp_imgs_PATH = get_new_temp_imgs_folder(0)
    OUTPUT = go.Output(DB, args.output_path, args.flux_balance_analysis, args.knockouts)
    if args.inchidb:
        verbose_print(args.verbose, 'STATUS: {} tanimoto threshold being used'.format(float(args.tanimoto_threshold)*100))
        cytosol_compartmentID = get_compartmentID_from_db(DB, 'cytosol')
        extracell_compartmentID = get_compartmentID_from_db(DB, 'extracellular')
        SIM = ss.TanimotoStructureSimilarity(R.targets, DB.get_all_compounds(),
                                             cytosol_compartmentID, extracell_compartmentID,
                                             args.tanimoto_threshold)
        return(SIM.finaltargets, R.ignorerxns, OUTPUT, temp_imgs_PATH)
    else:
        return(R.targets, R.ignorerxns, OUTPUT, temp_imgs_PATH)
    DB.conn.close()

def retrieve_database_info(args):
    '''
    Generates database or uses previously generated database.
    Can also add metacyc database to a Kbase metabolic database
    '''
    def get_compartment(DB):
        compartmentID_array = DB.get_compartment('cytosol')
        if compartmentID_array is None:
            compartmentID = 'c0'
        else:
            if compartmentID_array[0] == '':
                compartmentID = 'c0'
            else:
                compartmentID = compartmentID_array[0]
        return (compartmentID)
 
    if args.generate_database:
        '''
        Generate a database
        '''
        init_db.Createdb(args.generate_database, args.inchidb)
        DB = Q.Connector(args.generate_database)
        if args.kbase:
            '''
            Add kbase daatabase
            '''
            bkdb.BuildKbase(args.kbase_dump_directory, PATH+'/Database/KbasetoKEGGCPD.txt',
                            PATH+'/Database/KbasetoKEGGRXN.txt', args.inchidb,
                            args.generate_database, args.kbase_reaction_type)
        if args.metacyc:
            '''
            Add metacyc daatabase
            '''
            bmcdb.Translate(args.generate_database, DB, args.metacyc_addition,args.translation_file,
                            args.inchidb, args.metacyc_reaction_type)
        if args.kegg and (args.kbase or args.metacyc):
            '''
            Add kegg daatabase
            '''
            BKD = bkeggdb.CompileKEGGIntoDB(args.generate_database, args.kegg_organism_type,
                                            args.inchidb, args.processors, args.kegg_number_of_organisms,
                                            args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, True)
            DB = BKD.DB

        elif args.kegg and not args.kbase and not args.metacyc:
            '''
            Add only kbase daatabase
            '''
            print ('STATUS: Only KEGG')
            BKD = bkeggdb.CompileKEGGIntoDB(args.generate_database, args.kegg_organism_type,
                                            args.inchidb, args.processors,
                                            args.kegg_number_of_organisms, args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, False)
            DB = BKD.DB

        if args.SPRESI:
            '''
            Translate synthetic rdf file to database
            '''
            compartmentID = get_compartment(DB)
            bspresidb.RDF_Reader(args.spresi_dump_directory,
                                      args.generate_database,
                                      args.spresi_reaction_type,
                                      compartmentID, args.processors)
            DB = Q.Connector(args.generate_database)

        if args.mine:
            bminedb.BuildMINEdb(args.mine_dump_directory, args.generate_database,
                                args.inchidb, args.mine_reaction_type)

        if args.atlas:
            batlasdb.build_atlas(args.atlas_dump_directory, args.generate_database, args.inchidb,
                                 args.processors, args.atlas_reaction_type)
        database = args.generate_database
        rdc.OverlappingCpdIDs(database)

    elif args.database:
        '''
        Use and existing database
        '''
        DB = Q.Connector(args.database)
        compartmentID = get_compartment(DB)

        if args.kbase:
            '''
            Add kbase daatabase
            '''
            bkdb.BuildKbase(args.kbase_dump_directory, PATH+'/Database/KbasetoKEGGCPD.txt', PATH+'/KbasetoKEGGRXN.txt',
                             args.inchidb, args.database, args.rxntype)

        if args.metacyc:
            '''
            Add metacyc daatabase
            '''     
            bmcdb.Translate(args.database, DB, args.metacyc_addition, args.translation_file,
                            args.inchidb, args.metacyc_reaction_type)

        if args.kegg:
            '''
            Add KEGG daatabase
            '''
            BKD = bkeggdb.CompileKEGGIntoDB(args.database, args.kegg_organism_type,
                                            args.inchidb, args.processors, args.kegg_number_of_organisms,
                                            args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, True)
            DB = BKD.DB

        if args.SPRESI:
            '''
            Add a synthetic rdf file to the database
            '''
            bspresidb.RDF_Reader(args.spresi_dump_directory,
                                      args.database,
                                      args.spresi_reaction_type,
                                      compartmentID, args.processors)
            DB = Q.Connector(args.database)

        if args.mine:
            bminedb.BuildMINEdb(args.mine_dump_directory, args.database,
                                args.inchidb, args.mine_reaction_type)

        if args.atlas:
            batlasdb.build_atlas(args.atlas_dump_directory, args.database, args.inchidb,
                                 args.processors, args.atlas_reaction_type)
        if args.inchidb and (args.kbase or args.metacyc or args.kegg or args.SPRESI or args.mine or args.atlas):
            '''only run to remove overlapping ids if a database was added''' 
            rdc.OverlappingCpdIDs(database)

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
                                   reverse_constraints=False)
        with open(args.generate_database_constraints, 'wb') as fout1:
            cPickle.dump(LP.A, fout1)
            cPickle.dump(LP.allcpds, fout1)

    elif args.database_constraints:
        with open(args.database_constraints, 'rb') as fin1:
            A = cPickle.load(fin1)
            allcompounds4matrix = cPickle.load(fin1)
        LP = co.ConstructInitialLP(allrxns, allcompounds4matrix, DB,
                                   ignore_reactions, A,
                                   reverse_constraints=False)
    DB.conn.close()
    return LP

def construct_and_run_integerprogram(args, targets, output, database):
    '''
    Constructs ILP and solves it identifying shortest path to the target
    '''
    DB = Q.Connector(database)
    IP = ip_pulp.IntergerProgram(DB, args.limit_reactions,
                                 args.limit_cycles, args.k_number_of_paths, args.cycles,
                                 args.solver_time_limit, output)
    if args.flux_balance_analysis:
        active_metabolism = retrieve_active_FBA_metabolism(targets, DB, args, output)
    else:
        active_metabolism = {}
    if args.predict_toxicity:
        print ('WARNING: currently will work only with host organisms E. Coli')
        toxicity_train = tt.TrainToxModel()
        return (IP, active_metabolism, toxicity_train)
    else:
        return (IP, active_metabolism, None)

def _specific_target(target_id):
    '''Determines if there was a specified organism'''
    if target_id in ['', 'NA', 'N/A']:
        return False
    else:
        return True

def retrieve_shortestpath(target_info, IP, LP, database, args, output, active_metabolism, toxicity_train, temp_imgs_PATH):
    '''Retrieve the shortest path for target organism'''
    DB = Q.Connector(database)
    verbose_print(args.verbose, "STATUS: getting path for {}".format(target_info))
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
            optimal_pathways = IP.run_glpk(LP, incpds_active, inrxns, target_info[0],
                                           multiplesolutions=args.multiple_solutions)
            if optimal_pathways:                    
                uniq_externalrxns = []
                for path in optimal_pathways:
                    path_org = []
                    for rxn in path:
                        rxn = re.sub('_F$', '', rxn)
                        rxn = re.sub('_R$', '', rxn)
                        path_org.append(rxn)
                    uniq_externalrxns.append(list(set(path_org) - set(inrxns)))

                ex_info = ei.Extract_Information(optimal_pathways, incpds_active, DB)
                if toxicity_train:
                    tox_excpd = toxicity_train.predict_toxicity(ex_info.temp_exmets)
                    output.output_shortest_paths(target_info, ex_info.temp_rxns, tox_excpd)
                else:
                    output.output_shortest_paths(target_info, ex_info.temp_rxns)
                if args.flux_balance_analysis:
                    opt_fba = run_flux_balance_analysis(target_info, ex_info, incpds,
                                                        incpds_active, inrxns_active,
                                                        args.media_for_FBA, args.knockouts,
                                                        output, DB)
                    R = rf.ReactionFiles(args.output_path, DB, ex_info.temp_rxns,
                                     target_info[0], target_info[2], incpds_active, args.figures)
                    output.output_raw_solutions(target_info[0], target_info[2], R.ordered_paths, ex_info.temp_rxns, incpds_active)
                    if args.figures:
                        G = spgd.GraphDot(DB, args.output_path, incpds_active, inrxns,
                                          temp_imgs_PATH, opt_fba.fbasol.x_dict)
                        G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)

                elif args.figures and not args.flux_balance_analysis:
                    G = spgd.GraphDot(DB, args.output_path, incpds, inrxns, temp_imgs_PATH)
                    G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)
            else:
                output.output_shortest_paths(target_info, [])
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
        comparisonKOresults = crko.CompareKO(target_info[0], opt_fba.compounds_dict, opt_fba.fbasol, opt_fba.KOsolutions,
                                             ex_info.temp_rxns, DB)
        output.output_FBA_KOs(target_info, opt_fba.fbasol, opt_fba.compounds_dict, comparisonKOresults, ex_info.temp_rxns)
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
    targets, ignore_reactions, output, temp_imgs_PATH = read_in_and_generate_output_files(args, database)
    LP = retrieve_constraints(args, all_db_reactions, all_db_compounds, ignore_reactions, database)
    IP, active_metabolism, toxicity_train = construct_and_run_integerprogram(args, targets, output, database)

    def start_processes_new(index, target):
        verbose_print(args.verbose, 'STATUS: Inititate new process for target {}'.format(target))
        p = Process(target=retrieve_shortestpath, args=(target, IP, LP, database, args, output,
                                                        active_metabolism, toxicity_train, temp_imgs_PATH))
        p.start()
        return (p)

    def start_processes(targets, processors):
        verbose_print(args.verbose, 'STATUS: Initiating solving of solutions for initial of {} targets'.format(processors))
        processes = []
        for i in range(0, int(processors)):
            try:
                processes.append(Process(target=retrieve_shortestpath, args=(targets[i], IP, LP, database, args, output,
                                                                             active_metabolism, toxicity_train, temp_imgs_PATH)))
                index = i
            except IndexError:  
                pass
        for p in processes:
            p.start()
        sleep_time = float(processors)/float(10)
        while index < len(targets)-1:
            time.sleep(sleep_time)
            for p in processes:
                if not p.is_alive():
                    index+=1
                    if index < len(targets)-1:
                        processes.remove(p)
                        p_new = start_processes_new(index, targets[index])
                        processes.append(p_new)
                    else:
                        verbose_print(args.verbose, 'STATUS: All targets optimal solutions have been started')
                        break
        if index == len(targets)-1:
            verbose_print(args.verbose, 'STATUS: Solving last set of solutions for targets..waiting until all {} targets are complete'.format(len(processes)))
            last_target = True
            while last_target is True:
                time.sleep(sleep_time)
                for p in processes:
                    if not p.is_alive():
                        verbose_print(args.verbose, 'STATUS: Starting last target')
                        p_new = start_processes_new(index, targets[index])
                        last_target = False
                        break
            for p in processes:
                p.join()
        if index >= len(targets)-1:
            verbose_print(args.verbose, 'STATUS: Waiting for all target solutions')
            for p in processes:
                p.join()
 
    start_processes(targets, args.processors)
    if args.output_xlsx_format:
        output.convert_output_2_xlsx()

    '''Remove all temporary images'''
    if args.images == 'True':
        for filename in glob.glob(PATH+"/Visualization/compound*"):
            os.remove(filename)

    '''Removes all dot files if they exist'''
    for filename in glob.glob(args.output_path+"/solution_figures/*.dot"):
        os.remove(filename)

if __name__ == '__main__':
    main()
    exit()
