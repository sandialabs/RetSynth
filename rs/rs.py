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
import shutil
from timeit import default_timer as timer
from Parser import read_startcompounds as rtsc
from Parser import read_targets as rt
from Parser import generate_output as go
from Parser import structure_similarity as ss
from Visualization_chemdraw import reaction_files as rf
from ShortestPath import extractinfo as ei
from ShortestPath import constraints as co
from ShortestPath import integerprogram_pulp as ip_pulp
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import build_modelseed as bms
from Database import build_metacyc_db as bmcdb
from Database import build_user_rxns_db as burdb
from Database import build_ATLAS_db as batlasdb
from Database import build_MINE_db as bminedb
from Database import build_KEGG_db as bkeggdb
from Database import build_SPRESI_db as bspresidb
from Database import query as Q
from Database import remove_duplicate_cpds as rdc
from FBA import build_model as bm
from FBA import optimize_target as ot
from FBA import compare_results as cr
from FBA import retrieve_producable_mets as rpm
from FBA import compareKO_results as crko

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
    parser.add_argument('-t', '--targets', help='Input file containing target compounds and organisms \
                                                 (unless list of starting compounds are \
                                                 provided see --start_compounds option) \
                                                 (NOTE: file must be tab deliminated if users \
                                                 is specifying target and organism)',
                        required=True, type=str)
 
    parser.add_argument('-op', '--output_path', help='Destination for output files',
                        required=False, type=str, default=PATH)

    parser.add_argument('-oxf', '--output_xlsx_format', help='Convert output txt files to xlsx \
                                                              files (warning: do not do this if \
                                                              getting pathways for 100s of compounds)',
                        required=False, action="store_true")
    parser.add_argument('-stcpds', '--start_compounds', help='Insead of using an host organism FBA model \
                                                              as starting material \
                                                              the user can provide a \
                                                              a list of compounds that the shortest \
                                                              path should start from \
                                                              (use in place of target organism)',
                        required=False, type=str)
    ###NUMBER OF PROCESSOR AND VERBOSE OPTION
    parser.add_argument('-p', '--processors', help='Number of processors to use when \
                                                    solving for shortest path (default 4) \
                                                    should not exceed 100 \
                                                    (NOTE: many processors are only useful if getting \
                                                    pathways for many (~10) compounds)',
                        required=False, type=int, default=4)
    parser.add_argument('-v', '--verbose', help="Level of output", action="store_true")

    ###DATABASE OPTIONS###
    parser.add_argument('-gdb', '--generate_database', help='Generate RetSynth database', type=str)

    parser.add_argument('-gdbc', '--generate_database_constraints', help='Generate output constraint \
                                                                          file for entire database',
                        required=False, type=str)

    parser.add_argument('-db', '--database', help='Specify pre-built database to use',
                        required=False, type=str)

    parser.add_argument('-dbc', '--database_constraints', help='Utilize pre-constructed constraint \
                                                                file for entire database',
                        required=False, type=str)    



    parser.add_argument('--inchidb', help='Retrieve InChis and use them as compound \
                                           IDs in the metabolic database (NOTE: This extends the time \
                                           needed to build database but makes pathway identification \
                                           as compound IDs are more uniform)', 
                        action='store_true')

    ###INTEGRATING PATRIC FBA MODELS INTO RETSYNTH DB
    parser.add_argument('--patric_models', help='Set whether to build database with patric data (using patric, patric.org), \
                                                 (NOTE: To use patric one must have a patric account, patric username and password are required!!)',
                        required=False, action="store_true")

    parser.add_argument('-p_un', '--patric_username', help='Specify patric username',  required=False, type=str) 

    parser.add_argument('-p_pw', '--patric_password', help='Specify patric password',  required=False, type=str) 

    parser.add_argument('-patricrxntype', '--patric_reaction_type', help='Define type of reactions \
                                                                          that are being added to database \
                                                                          (options are bio (default) and chem)',
                        required=False, type=str, default='bio')

    parser.add_argument('--patric_media', help='Type of media to build FBA models on (Complete (default) \
                                                (NOTE: Media name must be an exact match to file name in /home/media folder)',
                        
                        required=False, type=str, default='Complete')

    parser.add_argument('--patricfile', help='List of patric genomes user wants to add into RetSynth database (Default is list includes Escherichia coli DH1, \
                                              Streptomyces venezuelae ATCC 10712 and Streptomyces venezuelae ATCC 10712 genomes) if user wants to use their \
                                              own specific set of models file must be a csv file with header followed by the the list of patric genome IDs and Names \
                                              (i.e genomeID,genomeName) (NOTE: Also patric file can be generated on patric website. To do this go to https://www.patricbrc.org and \
                                              go to Organisms banner and then click on All Bacteria.  Then click on the Genomes tab and select the genomes \
                                              you want to add to your database and then click the download select (csv) button',
                        required=False, type=str, default=PATH+'/Database/PATRIC_genome_complete_07152018.csv')

    parser.add_argument('--patric_sbml_output', help='Output FBA models generated by patric as SBML files \
                                                      (will be stored in the output_path folder)',
                        required=False, action="store_true")

    parser.add_argument('--patric_models_already_built', help='If user wants add patric FBA models to \
                                                               the RetSynth Database \
                                                               that were previously \
                                                               generated by RetSynth \
                                                               stored in the patric users /home/models/ folder \
                                                               then this option should be specified',
                        required=False, action="store_true")

    ###INTEGRATING KBASE FBA MODELS INTO RETSYNTH DB
    parser.add_argument('--kbase', help='Set whether to build database with Kbase data, \
                                         requires folder of fba models (SBML format) previously \
                                         downloaded from kbase website https://kbase.us/',  required=False,
                        action="store_true")
                    
    parser.add_argument('-k_dir', '--kbase_dump_directory', help='Path to folder \
                                                                  of SBML network files from kbase',
                        required=False, type=str)
    parser.add_argument('-kbaserxntype', '--kbase_reaction_type', help='Define type of reactions \
                                                             that are being added to database \
                                                            (options are bio (default) and chem)',
                        required=False, type=str, default='bio')
    ###INTEGRATING METACYC REACTIONS INTO RETSYNTH DB
    parser.add_argument('--metacyc', help='Set whether to build database with metacyc',  required=False,
                        action="store_true")    
    
    parser.add_argument('-mc', '--metacyc_addition', help='Add metacyc xml file to database \
                                                           can be downloaded from the metacyc website \
                                                           https://metacyc.org/ should be named \
                                                           (metabolic-reactions.xml)',
                        required=False, type=str)
    
    parser.add_argument('-mcrxntype', '--metacyc_reaction_type', help='Define type of reactions \
                                                             that are being added to database \
                                                            (options are bio (default) and chem)',
                        required=False, type=str, default='bio')

    ###INTEGRATING KEGG REACTIONS INTO RETSYNTH DB
    parser.add_argument('--kegg', help='Set whether to build database with Kegg data',  required=False,
                        action="store_true")

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

    ###INTEGRATING ATlAS REACTIONS INTO RETSYNTH DB
    parser.add_argument('--atlas', help='Set whether to build database with ATLAS of biochemistry data lcsb-databases.epfl.ch/atlas/ \
                                         (NOTE: user will have to submit a request through website to get files of ATLAS information)',
                        required=False, action="store_true")

    parser.add_argument('-a_dir', '--atlas_dump_directory', help='Path to folder \
                        of ATLAS files (.csv)', required=False, type=str)

    parser.add_argument('-atlasrxntype', '--atlas_reaction_type', help='Define type of reactions \
                                                                    from atlas files that are being \
                                                                    added to database (options \
                                                                    are bio (default) and chem)',
                        required=False, type=str, default='bio')

    ###INTEGRATING MINE REACTIONS INTO RETSYNTH DB
    parser.add_argument('--mine', help='Set whether to build database with Metabolic In Silico Network Expansion (MINE) database \
                                         http://minedatabase.mcs.anl.gov/#/home data (NOTE: This option reads in msp files which can \
                                         be found on the website in the Download section)',
                        required=False, action="store_true")
    
    parser.add_argument('-m_dir', '--mine_dump_directory', help='Path to downloaded folder of .msp files from MINE database',
                        required=False, type=str)

    parser.add_argument('-minerxntype', '--mine_reaction_type', help='Define type of reactions \
                                                                      from atlas files that are being \
                                                                      added to database (options \
                                                                      are bio (default) and chem)',
                        required=False, type=str, default='bio')

    ###INTEGRATING SPRESI REACTIONS INTO RETSYNTH DB
    parser.add_argument('--SPRESI', help='Set whether to build database with SPRESI information \
                                          (NOTE: to integrate spresi database into RetSynth database \
                                          the user must have purchased this database and have it stored \
                                          in rdf file format)',
                        required=False, action="store_true")
    
    parser.add_argument('-s_dir', '--spresi_dump_directory', help='Path to folder of \
                                                                   SPRESI files (.rdf)',
                         required=False, type=str)

    parser.add_argument('-spresirxntype', '--spresi_reaction_type', help='Define type of reactions \
                                                                          from rdf files that are being \
                                                                          added to database (options \
                                                                          are bio and chem (default))',
                        required=False, type=str, default='chem')

    ###INTEGRATING USER IDENTIFIED INDIVIDUAL REACTIONS INTO RETSYNTH DB
    parser.add_argument('-user_rxns', '--user_rxns_2_database', help='Add user defined reactions  \
                                                                      to the metabolic database \
                                                                      (tab delimnated text file)',
                        required=False, type=str)
    parser.add_argument('-user_rxns_type', '--user_rxns_2_database_type', help='Define type of reactions \
                                                                                that are being added to database \
                                                                                (options are bio (default) and chem)',
                        required=False, type=str, default='bio')

    ###FLUX BALANCE ANALYSIS OPTIONS###    
    parser.add_argument('-fba', '--flux_balance_analysis', help='Runs flux balance analysis using cobrapy \
                                                                 on model organism of interest',
                        required=False, action="store_true")

    parser.add_argument('-media', '--media_for_FBA', help='Define media for Flux Balance Analysis \
                                                           Complete media is the default, has to match \
                                                           --patric_media, or if using precompiled database \
                                                           the --patric_media option that was used to build \
                                                           the database',
                        required=False, default='Complete', type=str)

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

    parser.add_argument('-tmlim', '--solver_time_limit', help='time limit for solver to solve shortest \
                                                                path (note: only function with PULP python \
                                                                solver package)',
                        required=False, type=str, default=30)

    parser.add_argument('-evalrxns', '--evaluate_reactions', help='Defines which type of reactions \
                                                                   (bio, chem, or all (default)) \
                                                                   to be evaluated in identifying \
                                                                   shortest paths',
                        required=False, type=str, default='all')
    
    ###EXTRA OPTIONS FOR THE USER TO SPECIFY THE TO GET THE TYPES OF PATHWAYS THEY NEED###
    parser.add_argument('-k', '--k_number_of_paths', help='Specifies the level of shortest \
                                                           pathways wanted by the user',
                        required=False, type=int, default=0)

    parser.add_argument('-ms', '--multiple_solutions', help='Find all shortest paths (True \
                                                             (defualt) or False)',
                        required=False, type=str, default='True')

    parser.add_argument('-cy', '--cycles', help='Elimate cycles when finding shortest paths \
                                                 (True (default) or False)',
                        required=False, type=str, default='True')

    parser.add_argument('-run_tan_thresh', '--run_tanimoto_threshold', help='Tells program to run tanimoto threshold analysis \
                                                                             identify strucutrally similar compounds in database to targets',
                         required=False, action='store_true')

    parser.add_argument('-tan_thresh', '--tanimoto_threshold', help='Set the tanimoto threshold to \
                                                                     identify related compounds',
                         required=False, type=float, default=1)


    ###FIGURE OPTIONS###
    parser.add_argument('--figures_chemdraw', help='Generate pathway chemdraw figures', action='store_true')

    parser.add_argument('--figures_graphviz', help='Generate pathway figures using graphviz', action='store_true')

    parser.add_argument('--show_chemical_rxn_info', help='In chemdraw figures show all chemical reaction information',
                        action='store_true')
    parser.add_argument('--images', help='Set whether to use chemical images \
                                          or round nodes in output figures\
                                          (default True)', required=False,
                        type=str, default='True')
    parser.add_argument('--timer_output', help='Generates output file calculating time \
                                                it took to do various functions througout \
                                                the software (default False)', required=False,
                        action="store_true")
    return parser.parse_args()



def check_arguments(args):
    '''Checks and makes sure all required arguments are provided'''
    parser = argparse.ArgumentParser()
    if args.kbase and not args.kbase_dump_directory:
        parser.error('Requires use of --kbase_dump_directory')

    if args.knockouts and not args.flux_balance_analysis:
        parser.error('--knockouts option requires that \
                     --flux_balance_analysis option be also specified')
    if args.database and (not args.generate_database_constraints and not args.database_constraints):
        print ('WARNING: User specified specific database but not constraint file therefore default constraint file will be used but may not match user specified database')
 
    if args.metacyc and not args.metacyc_addition:
        parser.error('--metacyc requires use of parameters metacyc_addition')        

    if args.atlas and not args.atlas_dump_directory:
        parser.error('--atlas requires use of --atlas_dump_directory')

    if args.SPRESI and not args.spresi_dump_directory:
        parser.error('--SPRESI requires use of --spresi_dump_directory')

    if args.mine and not args.mine_dump_directory:
        parser.error('--mine requires use of --mine_dump_directory')

    if args.patric_models and not args.patric_password:
        parser.error('--patric_models requires options --patric_models, --patric_username, and --patric_password be specified')
    
    if args.patric_models and not args.patric_username:
        parser.error('--patric_models requires options --patric_models, --patric_username, and --patric_password be specified')
    
    if args.patric_username and not args.patric_password:
        parser.error('--patric_username requires options --patric_models, --patric_username, and --patric_password be specified')

    if not args.patric_models and args.patric_models_already_built:
        parser.error('--args.previously_built_patric_models requires options --patric_models, --patric_username, and --patric_password be specified')

    # if not args.patric_models and args.patricfile:
    #     parser.error('--args.patricfile requires options --patric_models, --patric_username, and --patric_password be specified')     

    if args.start_compounds and args.flux_balance_analysis:
        parser.error('Flux balance cannot be performed on a set of starting compounds\
                     would need to use an organisms metabolism to simulate flux')

    if not args.multiple_solutions and args.k_number_of_paths:
        parser.error('Cannot find k_number_of_paths correctly \
                     without finding all multiple_solutions')
    if not args.targets:
        parser.error('Requires an input file of target compounds')

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

def get_new_temp_imgs_folder(PATH, count):
    '''Check if folder to store images is already present if so new temp folder is made'''
    count+=1
    try:
        os.mkdir(PATH+'/temp_imgs_'+str(count))
        return PATH+'/temp_imgs_'+str(count)
    except OSError:
        PATH_NEW = get_new_temp_imgs_folder(PATH, count)
        return PATH_NEW

def read_in_and_generate_output_files(args, database):
    '''Read in target input file and generate output files'''
    DB = Q.Connector(database)
    R = rt.Readfile(args.targets, DB, args.inchidb)
    if not R.targets:
        raise ValueError('ERROR: No targets, try different compounds')
    temp_imgs_PATH = get_new_temp_imgs_folder(args.output_path, 0)
    OUTPUT = go.Output(DB, args.output_path, args.verbose, args.flux_balance_analysis, args.knockouts, args.timer_output)
    if args.run_tanimoto_threshold:
        verbose_print(args.verbose, 'STATUS: {} tanimoto threshold being used'.format(float(args.tanimoto_threshold)*100))
        cytosol_compartmentID = get_compartmentID_from_db(DB, 'cytosol')
        extracell_compartmentID = get_compartmentID_from_db(DB, 'extracellular')
        print (args.tanimoto_threshold)
        SIM = ss.TanimotoStructureSimilarity(R.targets, DB.get_all_compounds(),
                                             cytosol_compartmentID, extracell_compartmentID,
                                             args.verbose, args.tanimoto_threshold)
        OUTPUT.output_final_targets(SIM.finaltargets, args.tanimoto_threshold)
        return(SIM.finaltargets, R.ignorerxns, OUTPUT, temp_imgs_PATH)

    else:
        return(R.targets, R.ignorerxns, OUTPUT, temp_imgs_PATH)

def retrieve_database_info(args):
    '''
    Generates database or uses previously generated database.
    Can also add metacyc database to a Kbase metabolic database
    '''
    if args.generate_database:
        '''
        Generate a database
        '''
        init_db.Createdb(args.generate_database, args.inchidb)
        DB = Q.Connector(args.generate_database)
        if args.patric_models:
            bms.BuildModelSeed(username=args.patric_username, password=args.patric_password, rxntype=args.patric_reaction_type,
                               inchidb=args.inchidb, DBpath=args.generate_database, output_folder=args.output_path, media=args.patric_media, 
                               patricfile=args.patricfile, newdb=True, sbml_output=args.patric_sbml_output,
                               previously_built_patric_models=args.patric_models_already_built)
        if args.kbase:
            '''
            Add kbase daatabase
            '''
            bkdb.BuildKbase(args.kbase_dump_directory, PATH+'/Database/data/KbasetoKEGGCPD.txt',
                            PATH+'/Database/data/KbasetoKEGGRXN.txt', args.inchidb,
                            args.generate_database, args.kbase_reaction_type)
        if args.metacyc:
            '''
            Add metacyc daatabase
            '''
            bmcdb.Translate(args.generate_database, args.metacyc_addition,
                            args.inchidb, args.metacyc_reaction_type, args.verbose)
        if args.kegg and (args.patric_models or args.kbase or args.metacyc):
            '''
            Add kegg daatabase
            '''
            BKD = bkeggdb.CompileKEGGIntoDB(args.generate_database, args.kegg_organism_type,
                                            args.inchidb, args.processors, args.kegg_number_of_organisms,
                                            args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, True)
            DB = BKD.DB

        elif args.kegg and not args.kbase and not args.patric_models and not args.metacyc:
            '''
            Add only kbase daatabase
            '''
            print ('STATUS: Only KEGG')
            BKD = bkeggdb.CompileKEGGIntoDB(args.generate_database, args.kegg_organism_type,
                                            args.inchidb, args.processors,
                                            args.kegg_number_of_organisms, args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, False)
            DB = BKD.DB

        if args.user_rxns_2_database:
            burdb.AddUserRxns2DB(args.generate_database, args.user_rxns_2_database,
                                 model_id='UserAdded', rxntype=args.user_rxns_2_database_type)
        if args.SPRESI:
            '''
            Translate synthetic rdf file to database
            '''
            cytosol_compartmentID = get_compartmentID_from_db(DB, 'cytosol')
            bspresidb.RDF_Reader(args.spresi_dump_directory,
                                      args.generate_database,
                                      args.spresi_reaction_type,
                                      cytosol_compartmentID, args.processors)
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
        if args.patric_models:
            bms.BuildModelSeed(username=args.patric_username, password=args.patric_password, rxntype=args.patric_reaction_type,
                               inchidb=args.inchidb, DBpath=args.generate_database, output_folder=args.output_path, media=args.patric_media, 
                               patricfile=args.patricfile, newdb=False, sbml_output=args.patric_sbml_output, 
                               previously_built_patric_models=args.patric_models_already_built)
        if args.kbase:
            '''
            Add kbase daatabase
            '''
            bkdb.BuildKbase(args.kbase_dump_directory, PATH+'/Database/data/KbasetoKEGGCPD.txt',
                            PATH+'/Database/data/KbasetoKEGGRXN.txt',
                            args.inchidb, args.database, args.rxntype)

        if args.metacyc:
            '''
            Add metacyc daatabase
            '''     
            bmcdb.Translate(args.database, args.metacyc_addition,
                            args.inchidb, args.metacyc_reaction_type, args.verbose)

        if args.kegg and (args.patric_models or args.kbase or args.metacyc):
            '''
            Add kegg daatabase
            '''
            BKD = bkeggdb.CompileKEGGIntoDB(args.generate_database, args.kegg_organism_type,
                                            args.inchidb, args.processors, args.kegg_number_of_organisms,
                                            args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, True)
            DB = BKD.DB

        elif args.kegg and not args.kbase and not args.patric_models and not args.metacyc:
            '''
            Add only kbase daatabase
            '''
            print ('STATUS: Only KEGG')
            BKD = bkeggdb.CompileKEGGIntoDB(args.generate_database, args.kegg_organism_type,
                                            args.inchidb, args.processors,
                                            args.kegg_number_of_organisms, args.kegg_number_of_organism_pathways,
                                            args.kegg_reaction_type, False)
            DB = BKD.DB

        if args.user_rxns_2_database:
            burdb.AddUserRxns2DB(args.database, args.user_rxns_2_database,
                                 model_id='UserAdded', rxntype=args.user_rxns_2_database_type)
        if args.SPRESI:
            '''
            Add a synthetic rdf file to the database
            '''
            cytosol_compartmentID = get_compartmentID_from_db(DB, 'cytosol')
            bspresidb.RDF_Reader(args.spresi_dump_directory,
                                      args.database,
                                      args.spresi_reaction_type,
                                      cytosol_compartmentID, args.processors)
            DB = Q.Connector(args.database)

        if args.mine:
            bminedb.BuildMINEdb(args.mine_dump_directory, args.database,
                                args.inchidb, args.mine_reaction_type)

        if args.atlas:
            batlasdb.build_atlas(args.atlas_dump_directory, args.database, args.inchidb,
                                 args.processors, args.atlas_reaction_type)
        database = args.database
        if args.inchidb and (args.patric_models or args.kbase or args.metacyc or args.kegg or args.SPRESI or args.mine or args.atlas):
            '''only run to remove overlapping ids if a database was added''' 
            rdc.OverlappingCpdIDs(database)

    else:
        if args.media_for_FBA == 'Carbon-D-Glucose':
            print ('WARNING: No database specified using pre constructed database with media {}'.format(args.media_for_FBA))
            database = PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_GL_MC.db'
            DB = Q.Connector(database)
        elif args.media_for_FBA == 'Complete':
            print ('WARNING: No database specified using pre constructed database with media {}'.format(args.media_for_FBA))
            
            # database = PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC_SPRESI.db'
            # database = PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC.db'
            database = PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC_cas.db'
            
            DB = Q.Connector(database)

    allcpds = DB.get_all_compounds()
    if args.evaluate_reactions == 'all':
        allrxns = DB.get_all_reactions()
    elif args.evaluate_reactions == 'bio':
        allrxns = DB.get_reactions_based_on_type('bio')
    elif args.evaluate_reactions == 'chem':
        allrxns = DB.get_reactions_based_on_type('chem')
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
    else:
        def load_preconstructed_constraint_files(media, db, args):
            print ('WARNING: No database constraint file specified using pre constructed database constraint file for database with media {}'.format(args.media_for_FBA))
            if args.evaluate_reactions =='all' or args.evaluate_reactions =='bio':
                # with open(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC_SPRESI.constraints'.format(db), 'rb') as fin1:
                # with open(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC.constraints'.format(db), 'rb') as fin1:
                with open(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC_cas.constraints'.format(db), 'rb') as fin1:                
                    A = cPickle.load(fin1)
                    allcompounds4matrix = cPickle.load(fin1)
                return (A, allcompounds4matrix)
    
            elif args.evaluate_reactions =='chem':
                # print ('WARNING: At the moment pre constructed databases do not have chem reactions try bio reactions instead')
                # return (None, None)
                with open(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC_SPRESI_chem.constraints'.format(db), 'rb') as fin1:
                    A = cPickle.load(fin1)
                    allcompounds4matrix = cPickle.load(fin1)
                return (A, allcompounds4matrix)

        if args.media_for_FBA=='Carbon-D-Glucose':
            A, allcompounds4matrix = load_preconstructed_constraint_files(args.media_for_FBA, 'GL', args)
        elif args.media_for_FBA=='Complete':
            A, allcompounds4matrix = load_preconstructed_constraint_files(args.media_for_FBA, 'CP', args)
        if A and allcompounds4matrix:
            LP = co.ConstructInitialLP(allrxns, allcompounds4matrix, DB,
                                       ignore_reactions, A,
                                       reverse_constraints=False)
        else: 
            print ('ERROR: No identified pre constraint file...stopping run')
    return LP

def construct_and_run_integerprogram(args, targets, output, database):
    '''
    Constructs ILP and solves it identifying shortest path to the target
    '''
    DB = Q.Connector(database)
    if args.timer_output:
        IP = ip_pulp.IntergerProgram(DB, args.limit_reactions,
                                    args.limit_cycles, args.k_number_of_paths,
                                    args.cycles, args.verbose, args.solver_time_limit, output)
    else:
        IP = ip_pulp.IntergerProgram(DB, args.limit_reactions,
                                     args.limit_cycles, args.k_number_of_paths,
                                     args.cycles, args.verbose, args.solver_time_limit, args.timer_output)
    
    return (IP)

def _specific_target(target_id):
    '''Determines if there was a specified organism'''
    if target_id in ['', 'NA', 'N/A']:
        return False
    else:
        return True

def retrieve_shortestpath(target_info, IP, LP, database, args, output, temp_imgs_PATH):
    '''Retrieve the shortest path for target organism'''
    start = timer()
    DB = Q.Connector(database)
    verbose_print(args.verbose, "STATUS: getting path for {}".format(target_info))
    if args.images == 'False':
        _images = False
    else:
        _images = True
    if not _specific_target(target_info[2]) and not args.start_compounds:
        print ('WARNING: No organism given therefore target {} compound will be skipped ... '.format(target_info[0]))
    else:
        if args.start_compounds:
            incpds_active = rtsc.readfile_startcompounds(args.start_compounds)
            inrxns_active = []
        else:
            incpds_active = DB.get_compounds_in_model(target_info[2])
            inrxns_active = DB.get_reactions_in_model(target_info[2])

        if target_info[0] in incpds_active: #Check if compound exists in organism
            output.output_compound_natively_present_in_target_organism(target_info)
        else:
            optimal_pathways = IP.run_glpk(LP, incpds_active, inrxns_active, target_info[0],
                                           multiplesolutions=args.multiple_solutions)
            if optimal_pathways:                    
                uniq_externalrxns = []
                for path in optimal_pathways:
                    path_org = []
                    for rxn in path:
                        rxn = re.sub('_F$', '', rxn)
                        rxn = re.sub('_R$', '', rxn)
                        path_org.append(rxn)
                    uniq_externalrxns.append(list(set(path_org) - set(inrxns_active)))

                ex_info = ei.Extract_Information(optimal_pathways, incpds_active, inrxns_active, DB)
                output.output_shortest_paths(target_info, ex_info.temp_rxns)
        
                R = rf.ReactionFiles(args.output_path, DB, ex_info.temp_rxns,
                                 target_info[0], target_info[2], incpds_active, args.figures_graphviz)
                output.output_raw_solutions(target_info[0], target_info[2], R.ordered_paths,
                                            ex_info.temp_rxns, ex_info.temp_external, incpds_active)
                if args.flux_balance_analysis:
                    opt_fba = run_flux_balance_analysis(target_info, ex_info,
                                                        incpds_active, inrxns_active,
                                                        args.media_for_FBA, args.knockouts,
                                                        output, DB, args.verbose)

                    if args.figures_graphviz:
                        from Visualization_graphviz import SP_Graph_dot as spgd
                        G = spgd.GraphDot(DB, args.output_path, incpds_active, inrxns_active,
                                          temp_imgs_PATH, opt_fba.fbasol.fluxes)
                        G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)
                    if args.figures_chemdraw:
                        R.generate_cdxml_files(fba_fluxes=opt_fba.fbasol.fluxes, show_rxn_info=args.show_chemical_rxn_info)


                elif (args.figures_graphviz or args.figures_chemdraw) and not args.flux_balance_analysis:
                    if args.figures_graphviz:
                        from Visualization_graphviz import SP_Graph_dot as spgd
                        G = spgd.GraphDot(DB, args.output_path, incpds_active, inrxns_active, temp_imgs_PATH)
                        G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)
                    if args.figures_chemdraw:
                        R.generate_cdxml_files(show_rxn_info=args.show_chemical_rxn_info)
            else:
                output.output_shortest_paths(target_info, [])
                if args.flux_balance_analysis:
                    verbose_print(args.verbose, 'WARNING: No optimal path for %s in species %s therefore no flux balance will be performed' % (target_info[0], target_info[2]))
    end = timer()
    if args.timer_output:
       output.output_timer('Time to find all paths for {}\t{}\t{}\n'.format(target_info[0], (end-start), (end-start)/60))
    verbose_print(args.verbose, "Time to find all paths for "+str(target_info[0])+' '+str(end - start))

def run_flux_balance_analysis(target_info, ex_info, incpds_active,
                              inrxns, media, ko,
                              output, DB, verbose):
    '''
    Run flux balance analysis on target organism with added reactions
    necessary to produce target compound
    '''
    fba = bm.BuildModel(target_info[2], incpds_active, inrxns, DB, verbose, media)
    opt_fba = ot.OptimizeTarget(target_info[0], target_info[2], fba.model, ex_info.temp_rxns,
                                ex_info.temp_exmets, fba.compounds_dict, incpds_active,
                                inrxns, DB, verbose, ko)
    # print (opt_fba.fbasol.fluxes)
    comparisonresults = cr.Compare(target_info[0], fba.solution, opt_fba.fbasol,
                                   ex_info.temp_rxns, DB)
    output.output_FBA(target_info, fba.solution, opt_fba, comparisonresults, ex_info.temp_rxns)
    output.output_theoretical_yield(target_info[0], target_info[2], opt_fba.fbasol,
                                    opt_fba.compounds_dict)
    if ko:
        output.output_essential_reactions(target_info[0], target_info[2], opt_fba.essentialrxns)
        comparisonKOresults = crko.CompareKO(target_info[0], opt_fba.compounds_dict, opt_fba.fbasol,
                                             opt_fba.KOsolutions, ex_info.temp_rxns, DB)
        output.output_FBA_KOs(target_info, opt_fba.fbasol, opt_fba.compounds_dict, comparisonKOresults, ex_info.temp_rxns)
    return opt_fba


def main():
    '''Main class'''
    args = parse_arguments()
    check_arguments(args)
    all_db_compounds, all_db_reactions, database = retrieve_database_info(args)
    targets, ignore_reactions, output, temp_imgs_PATH = read_in_and_generate_output_files(args, database)
    LP = retrieve_constraints(args, all_db_reactions, all_db_compounds, ignore_reactions, database)
    IP = construct_and_run_integerprogram(args, targets, output, database)


    def start_processes_new(index, target):
        '''Start solving a new solution for a target compound'''
        verbose_print(args.verbose, 'STATUS: Inititate new process for target {}'.format(target))
        p = Process(target=retrieve_shortestpath, args=(target, IP, LP, database, args, output,
                                                        temp_imgs_PATH))
        p.start()
        return (p)

    def start_processes(targets, processors):
        '''Start solving for solutions '''
        verbose_print(args.verbose, 'STATUS: Initiating solving of solutions for initial of {} targets'.format(processors))
        processes = []
        for i in range(0, int(processors)):
            try:
                processes.append(Process(target=retrieve_shortestpath, args=(targets[i], IP, LP, database, args, output, temp_imgs_PATH)))
                index = i
            except IndexError:  
                pass
        for p in processes:
            p.start()
        sleep_time = float(processors)/float(10)
        while index <= len(targets)-1:
            time.sleep(sleep_time)
            for p in processes:
                if not p.is_alive():
                    index+=1
                    if index <= len(targets)-1:
                        processes.remove(p)
                        p_new = start_processes_new(index, targets[index])
                        processes.append(p_new)
                    else:
                        verbose_print(args.verbose, 'STATUS: All targets optimal solutions have been started')
                        break
        if index > len(targets)-1:
            verbose_print(args.verbose, 'STATUS: Waiting for all target solutions')
            for p in processes:
                p.join()

    start_processes(targets, args.processors)
    # for target in targets:
    #     retrieve_shortestpath(target, IP, LP, database, args, output,)

    if args.output_xlsx_format:
        output.convert_output_2_xlsx()

    '''Remove all temporary images'''
    shutil.rmtree(temp_imgs_PATH)

    '''Removes all dot files if they exist'''
    for filename in glob.glob(args.output_path+"/solution_figures/*.dot"):
        os.remove(filename)

if __name__ == '__main__':
    main()
