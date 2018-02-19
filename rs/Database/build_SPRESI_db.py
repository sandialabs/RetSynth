from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'RDF reader (specifically for SPRESI database)'

from multiprocessing import Pool
import sqlite3
import os
import re
import glob
import shutil
from sys import platform
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_inchi
elif platform == "win32" or platform == 'win64':
    raise ImportError('Cannot translate RDF file on windows machine')

PATH = os.path.dirname(os.path.abspath(__file__))

def RDF_Reader(file_directory, DBpath, rxntype, compartment, processors,
               temp_option=False, pressure_option=False, yield_option=False,
               time_option=False, catalyst_option=False, solvent_option=False):
    '''
    Adds data from RDF files into database (specifically works with spresi formated rdf files)
    '''
    options = [temp_option, pressure_option, yield_option, time_option,
               catalyst_option, solvent_option, compartment]
    rdf_files = glob.glob(os.path.join(file_directory, '*'))
    args = [(i, str(c), options) for c, i in enumerate(rdf_files)]

    try:
        os.mkdir(PATH+'/temp')
    except OSError:
        shutil.rmtree(PATH+'/temp')
        os.mkdir(PATH+'/temp')
    if processors > 1:
        pool = Pool(processes=processors)
        pool.map(open_file, args)
        pool.close()
        pool.join()
    else:
        for arg in args:
            open_file(arg)
    add_info_2_database(DBpath, rxntype, compartment)

def open_file(args):
    '''
    Opens in an RDF file
    '''
    file_name = args[0]
    filenumber = args[1]
    options = args[2]
    RDF_dict = {}
    print ('STATUS: Loading/Parsing in RDF file ...'.format(file_name))
    output_file = open(PATH+'/temp/output_'+str(filenumber)+'.txt', 'w')
    open_file = open(file_name)
    for count, line in enumerate(open_file):
        line = line.strip('\n\r') #Need to keep front spaces to get smiles from mol structure
        if count == 0:
            if line != '$RDFILE 1':
                print ('WARNING: File did not have appropriate start title,\
                        may not be read correctly')
            else:
                pass
        elif count == 1:
            pass
        else:
            if line.startswith('$RFMT'):
                if len(RDF_dict.keys()) > 0:
                    parse_file(output_file, RDF_dict, filenumber, file_name, options)
                    RDF_dict = {}
                    array = line.split()
                    RDF_dict[array[2]] = []
                else:
                    array = line.split()
                    RDF_dict[array[2]] = []
            else:
                RDF_dict[array[2]].append(line)
    parse_file(output_file, RDF_dict, filenumber, file_name, options)

def get_mol_structure(match, item, type_bool, result_dict, count_item=0, GET_RXN=True):
    '''
    retrieves rxn, solvents and catalyst information
    '''
    if GET_RXN:
        if item.startswith('$MOL'):
            count_item += 1
            result_dict[count_item] = []
            type_bool = True
    else:
        if item.endswith('MOLSTRUCTURE') and match is not None:
            count_item += 1
            result_dict[count_item] = []
            type_bool = True
    if item.startswith('M  END'):
        type_bool = False
    if type_bool:
        if item.startswith('$'):
            pass
        else:
            result_dict[count_item].append(item)

    return(result_dict, type_bool, count_item)

def get_data(match, datatype, item, result_array, type_bool, DATA_TYPE=None):
    '''
    retrieves a variety of other information for rxn (yield, reference etc...)
    '''
    if item.endswith(datatype):
        type_bool = True
    if item.startswith('$DTYPE') and match is None:
        type_bool = False
    if type_bool is True and match is None:
        item = re.sub('\$DATUM ', '', item)
        result_array.append(item)
    return(result_array, type_bool)

def get_complex_reference_details(item, reference_parameters, full_citation_string,
                                  reference_details_array, reference_details_bool):
    '''
    retrieves details for references with many parameters
    '''
    for reference_parameter in reference_parameters:
        data_type = reference_parameter
        reference_match = re.search(data_type, item)
        reference_details_array, reference_details_bool = get_data(reference_match, data_type,
                                                                   item, reference_details_array,
                                                                   reference_details_bool)
        if reference_details_array:
            full_citation_string += ''.join(reference_details_array)

    return full_citation_string

def process_rxntext(rxntext_array):
    '''
    processes string containing reaction information
    '''
    temp = ''
    time = ''
    pressure = ''

    rxntext_string = ''.join(rxntext_array)
    match_stage = re.search('Stage', rxntext_string)
    if match_stage:
        time_array = []
        temp_array = []
        pressure_array = []

        stages = rxntext_string.split('|')
        for stage in stages:
            stage = re.sub('Stage \d+\: ', '', stage)
            values = stage.split(', ')
            count_pressure = 0
            count_temp = 0
            count_time = 0
            for value in values:
                if value.endswith('min') or value.endswith('h'):
                    count_time += 1
                    time_array.append(value)

                if value.endswith('degree') or value.endswith('k'):
                    count_temp += 1
                    temp_array.append(value)

                if value.endswith('atm'):
                    count_pressure += 1
                    pressure_array.append(value)
            if not count_time:
                time_array.append('None')
            if not count_temp:
                temp_array.append('None')
            if not count_pressure:
                pressure_array.append('None')

            time = ','.join(time_array)
            temp = ','.join(temp_array)
            pressure = ','.join(pressure_array)

    else:
        values = rxntext_string.split(', ')
        for value in values:
            if value.endswith('min') or value.endswith('h'):
                time = value
            if value.endswith('degree') or value.endswith('k'):
                temp = value
            if value.endswith('atm'):
                pressure = value
    return(temp, time, pressure)

def parse_file(output_file, RDF_dict, filenumber, file_name, options):
    '''
    Parses elements of an RDFile outputs them to new next file
    '''
    compartment = options[6]
    for key, items in RDF_dict.iteritems():
        key = 'rxn'+key+'_s'
        components = items[4].split()
        mol_bool = False
        solvent_bool = False
        catalyst_bool = False
        rxntext_bool = False
        yield_bool = False
        reference_type_bool = False
        reference_details_bool = False
        reference_t_details_bool = False
        reference_a_details_bool = False
        reference_j_details_bool = False
        reference_co_details_bool = False
        reference_n_details_bool = False
        reference_cl_details_bool = False
        reference_o_details_bool = False

        rxntext_array = []
        yield_array = []
        reference_type_array = []
        reference_t_details_array = []
        reference_a_details_array = []
        reference_j_details_array = []
        reference_co_details_array = []
        reference_n_details_array = []
        reference_cl_details_array = []
        reference_o_details_array = []
        reference_details_array = []
        mol = {}
        solvents = {}
        catalysts = {}
        full_citation_string = ''

        count_mol = 0
        count_solvents = 0
        count_catalysts = 0

        for item in items:
            #Get reaction reactants and products
            mol, mol_bool, count_mol = get_mol_structure(None, item, mol_bool, mol,
                                                         count_mol, GET_RXN=True)

            #Get reaction catalyst
            match_catalyst = re.search('CATALYST\(', item)
            catalysts, catalyst_bool, count_catalysts = get_mol_structure(match_catalyst, item,
                                                                          catalyst_bool,
                                                                          catalysts,
                                                                          count_catalysts,
                                                                          GET_RXN=False)

            #Get reaction solvents
            match_solvent = re.search('SOLVENT\(', item)
            solvents, solvent_bool, count_solvents = get_mol_structure(match_solvent, item,
                                                                       solvent_bool,
                                                                       solvents,
                                                                       count_solvents,
                                                                       GET_RXN=False)

            #Get reaction text information which contains info such as temp,
            #pressure, and time needed to perform reaction
            rxntext_match = re.search('RXNTEXT$', item)
            data_type = 'RXNTEXT'
            rxntext_array, rxntext_bool = get_data(rxntext_match, data_type, item,
                                                   rxntext_array, rxntext_bool)
            #Get Yield information
            yield_match = re.search('YIELD$', item)
            data_type = 'YIELD'
            yield_array, yield_bool = get_data(yield_match, data_type, item,
                                               yield_array, yield_bool)
            yield_string = ''.join(yield_array)

            #Get type of reference reaction was retrieved from
            reference_match = re.search('LITREF:TYPE$', item)
            data_type = 'LITREF:TYPE'
            reference_type_array, reference_type_bool = get_data(reference_match, data_type,
                                                                 item, reference_type_array,
                                                                 reference_type_bool)
            reference_type_string = ''.join(reference_type_array)

            #Depending on what type of reference was identified extract information
            if reference_type_string == 'JOURNAL':
                data_type = 'FULL_CITATION'
                reference_match = re.search(data_type, item)
                reference_details_array, reference_details_bool = get_data(reference_match,
                                                                           data_type,
                                                                           item,
                                                                           reference_details_array,
                                                                           reference_details_bool)
                full_citation_string = ''.join(reference_details_array)
                full_citation_string = 'JOURNAL ARTICLE '+full_citation_string

            if reference_type_string == 'PATENT':
                data_type = 'TITLE'
                reference_match = re.search(data_type, item)
                reference_t_details_array, reference_t_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_t_details_array,
                                                                               reference_t_details_bool)
                data_type = 'AUTHOR'
                reference_match = re.search(data_type, item)
                reference_a_details_array, reference_a_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_a_details_array,
                                                                               reference_a_details_bool)
                data_type = 'JOURNAL_YEAR'
                reference_match = re.search(data_type, item)
                reference_j_details_array, reference_j_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_j_details_array,
                                                                               reference_j_details_bool)
                data_type = 'NUMBER'
                reference_match = re.search(data_type, item)
                reference_n_details_array, reference_n_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_n_details_array,
                                                                               reference_n_details_bool)
                data_type = 'COUNTRY'
                reference_match = re.search(data_type, item)
                reference_co_details_array, reference_co_details_bool = get_data(reference_match,
                                                                                 data_type,
                                                                                 item,
                                                                                 reference_co_details_array,
                                                                                 reference_co_details_bool)

                data_type = 'CLASS'
                reference_match = re.search(data_type, item)
                reference_cl_details_array, reference_cl_details_bool = get_data(reference_match,
                                                                                 data_type,
                                                                                 item,
                                                                                 reference_cl_details_array,
                                                                                 reference_cl_details_bool)
                data_type = 'OWNER'
                reference_match = re.search(data_type, item)
                reference_o_details_array, reference_o_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_o_details_array,
                                                                               reference_o_details_bool)

            if (reference_type_string == 'COLLECTIVE VOLUME' or reference_type_string == 'DEPOTARTIKEL'
                    or reference_type_string == 'DISSERTATION' or reference_type_string == 'BUCH'):
                reference_parameters = ['TITLE', 'AUTHOR', 'JOURNAL_YEAR']
                data_type = 'TITLE'
                reference_match = re.search(data_type, item)
                reference_t_details_array, reference_t_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_t_details_array,
                                                                               reference_t_details_bool)
                data_type = 'AUTHOR'
                reference_match = re.search(data_type, item)
                reference_a_details_array, reference_a_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_a_details_array,
                                                                               reference_a_details_bool)
                data_type = 'JOURNAL_YEAR'
                reference_match = re.search(data_type, item)
                reference_j_details_array, reference_j_details_bool = get_data(reference_match,
                                                                               data_type,
                                                                               item,
                                                                               reference_j_details_array,
                                                                               reference_j_details_bool)

        #Print warning if reference type was not identified by code
        if reference_type_string == 'PATENT':
            full_citation_string = 'PATENT '+', '.join([''.join(reference_t_details_array),
            	                                           ''.join(reference_a_details_array),
                                                        ''.join(reference_j_details_array),
                                                        'Patent number-'+''.join(reference_n_details_array),
                                                        ''.join(reference_co_details_array),
                                                        'Patent Class-'+''.join(reference_cl_details_array),
                                                        'Patent Owner-'+''.join(reference_o_details_array)])

        elif (reference_type_string == 'COLLECTIVE VOLUME' or reference_type_string == 'DEPOTARTIKEL'
              or reference_type_string == 'DISSERTATION' or reference_type_string == 'BUCH'):
            full_citation_string = reference_type_string+' '+', '.join([''.join(reference_t_details_array),
                                                                        ''.join(reference_a_details_array),
                                                                        ''.join(reference_j_details_array)])
        elif not full_citation_string:
            print ('WARNING: Reaction has reference that is not patent or journal article: {} {}'.format(key, reference_type_string))

        #Determine whether to include reaction based on specified options
        temp, time, pressure = process_rxntext(rxntext_array)
        len_temp = len(temp.split(','))
        len_time = len(time.split(','))
        len_pressure = len(pressure.split(','))

        if (options[0] is True and temp == '') or (options[0] is True and len_temp == temp.count('None')):
            break

        if (options[1] is True and pressure == '') or (options[1] is True and len_pressure == pressure.count('None')):
            break

        if options[2] is True and yield_string == '':
            break

        if (options[3] is True and time == '') or (options[3] is True and len_time == time.count('None')):
            break

        if options[4] is True and len(catalysts) == 0:
            break
        if options[5] is True and len(solvents) == 0:
            break

        mol_smiles = generate_mol_file(mol, filenumber)
        catalysts_smiles = generate_mol_file(catalysts, filenumber)
        solvents_smiles = generate_mol_file(solvents, filenumber)

        key_error = False
        compound_error = False
        r = 0
        p = 0
        reactants = []
        products = []
        catalysts = []
        solvents = []
        for i in range(int(components[0])):
            r = i+1
            try:
                if mol_smiles[r] is not False:
                    reactants.append(str(mol_smiles[r])+'_'+compartment+'|---|'+'None')
                else:
                    compound_error = True
            except KeyError:
                key_error = True
                print ('WARNING: Reactant count off for reaction {} in {}'.format(key, file_name))

        for i in range(r, (int(components[1])+r)):
            p = i+1
            try:
                if mol_smiles[p] is not False:
                    products.append(str(mol_smiles[p])+'_'+compartment+'|---|'+'None')
                else:
                    compound_error = True
            except KeyError:
                key_error = True
                print ('WARNING: Product count off for reaction {} in {}'.format(key, file_name))

        for c in catalysts_smiles:
            try:
                if catalysts_smiles[c] is not False:
                    catalysts.append(str(catalysts_smiles[c])+'_'+compartment+'|---|'+'None')
                else:
                    compound_error = True
            except KeyError:
                key_error = True
                print ('WARNING: Catalyst count off for reaction {} in {}'.format(key, file_name))

        for s in solvents_smiles:
            try:
                if solvents_smiles[s] is not False:
                    solvents.append(str(solvents_smiles[s])+'_'+compartment+'|---|'+'None')
                else:
                    compound_error = True
            except KeyError:
                key_error = True
                print ('WARNING: Solvent count off for reaction {} in {}'.format(key, file_name))

        if key_error is False and compound_error is False:
            if temp == '':
                temp = None
            if pressure == '':
                pressure = None
            if time == '':
                time = None
            if yield_string == '':
                yield_string = None
            if full_citation_string == '':
                full_citation_string = None
            output_file.write('\t'.join([str(key), str("|----|".join(reactants)),
                                         str("|----|".join(products)),
                                         str("|----|".join(catalysts)),
                                         str("|----|".join(solvents)), str(temp), str(pressure),
                                         str(time), str(yield_string), str(full_citation_string)])+'\n')
        else:
            if compound_error is True:
                print ('{} skipped because of file issue'.format(key))
            if key_error is True:
                print ('{} skipped because of file issue'.format(key))

def generate_mol_file(compounds, filenumber):
    '''
    Generates mole file and then reads it in using the indigo API to get the smile
    '''
    IN = indigo.Indigo()
    INCHI = indigo_inchi.IndigoInchi(IN)
    compound_dict = {}
    for key, values in compounds.iteritems():
        output_mol_file = open(PATH+'/mol_output_'+str(filenumber)+'.mol', 'w')
        for value in values:
            output_mol_file.write(value+'\n')
        output_mol_file.close()
        try:
            mol = IN.loadMoleculeFromFile(PATH+'/mol_output_'+str(filenumber)+'.mol')
            try:
                inchi_value = INCHI.getInchi(mol)
                compound_dict[key] = inchi_value
            except indigo.IndigoException:
                try:
                    smile = mol.smiles()
                    compound_dict[key] = smile
                except indigo.IndigoException:
                    compound_dict[key] = False
            os.remove(PATH+'/mol_output_'+str(filenumber)+'.mol')
        except indigo.IndigoException:
            compound_dict[key] = False
    return compound_dict

def add_info_2_database(DBpath, rxntype, compartment):
    '''
    Adds data from text files too database
    '''
    conn = sqlite3.connect(DBpath, check_same_thread=False)
    conn.text_factory = str
    cnx = conn.cursor()
    cnx.execute("PRAGMA synchronous = OFF")
    cnx.execute("PRAGMA journal_mode = OFF")
    try:
        cnx.execute('''CREATE table reaction_catalysts
                    (reaction_ID text,catalysts_ID text,name text)''')
        cnx.execute('''CREATE INDEX reactioncatalysts_ind ON
                        reaction_catalysts(reaction_ID)''')
        cnx.execute('''CREATE table reaction_solvents
                    (reaction_ID text,solvents_ID text,name text)''')
        cnx.execute('''CREATE INDEX reactionsolvents_ind ON
                        reaction_solvents(reaction_ID)''')
        cnx.execute('''CREATE table reaction_spresi_info
                    (reaction_ID text, temperature text, pressure text,
                    total_time text, yield text, reference text)''')
        cnx.execute('''CREATE INDEX reaction_spresi_info_ind ON
                        reaction_spresi_info(reaction_ID)''')
    except sqlite3.OperationalError:
        pass
    Q = cnx.execute("SELECT ID FROM model WHERE ID = ?", ('SR1',))
    result = Q.fetchone()
    if result is None:
        Q = cnx.execute('''SELECT DISTINCT cluster_num FROM cluster''')
        hits = Q.fetchall()
        uniq_clusters = [i[0] for i in hits]
        cnx.execute("INSERT INTO cluster VALUES (?,?)", (len(uniq_clusters)+1, 'SR1'))
        cnx.execute("INSERT INTO model VALUES (?,?)", ('SR1', 'Synthetic_reactions'))

    text_files = glob.glob(os.path.join(PATH+'/temp', '*'))
    for text_file in text_files:
        print ('STATUS: Adding {} to database'.format(text_file))
        add_individual_file_info(text_file, cnx, conn, rxntype, compartment)

    cnx.execute("""DELETE FROM model_compound WHERE rowid NOT IN
                   (SELECT min(rowid) from model_compound group by cpd_ID,model_ID)""")
    cnx.execute("""DELETE FROM compound WHERE rowid NOT IN
                   (SELECT min(rowid) from compound group by ID)""")
    conn.commit()
    print ('STATUS: Reindexing indicies ...')
    cnx.execute('''REINDEX reactioncompound_ind1''')
    cnx.execute('''REINDEX reactioncompound_ind2''')
    cnx.execute('''REINDEX modelreaction_ind1''')
    cnx.execute('''REINDEX modelreaction_ind2''')
    cnx.execute('''REINDEX modelcompound_ind1''')
    cnx.execute('''REINDEX modelcompound_ind2''')
    cnx.execute('''REINDEX model_ind''')
    cnx.execute('''REINDEX reaction_ind''')
    cnx.execute('''REINDEX compound_ind''')
    cnx.execute('''REINDEX reactiongene_ind''')
    cnx.execute('''REINDEX reactionprotein_ind''')
    cnx.execute('''REINDEX cluster_ind''')
    cnx.execute('''REINDEX reaction_reversibility_ind''')
    cnx.execute('''REINDEX reactioncatalysts_ind''')
    cnx.execute('''REINDEX reactionsolvents_ind''')
    cnx.execute('''REINDEX reaction_spresi_info_ind''')

    try:
        cnx.execute('''REINDEX original_db_cpdIDs_ind''')
    except sqlite3.OperationalError:
        pass
    conn.commit()

    shutil.rmtree(PATH+'/temp')

def add_individual_file_info(text_file, cnx, conn, rxntype, compartment):
    '''
    Specifically adds individual file info from to database
    '''
    modelcompounds = [] #compoundID, #org
    allcompounds = [] #compoundID, name, compartment
    with open(text_file) as open_file:
        for line in open_file:
            line = line.strip('\n\r')
            larray = line.split('\t')
            cnx.execute("INSERT INTO model_reaction VALUES (?,?,?)", (larray[0], 'SR1', 'false'))
            cnx.execute("INSERT INTO reaction VALUES (?,?, ?, ?)", (larray[0], 'None', 'None', rxntype))
            cnx.execute("INSERT INTO reaction_reversibility VALUES (?,?)", (larray[0], 'false'))
            cnx.execute("INSERT INTO reaction_gene VALUES (?,?,?)", (larray[0], 'SR1', 'None'))
            cnx.execute("INSERT INTO reaction_protein VALUES (?,?,?)", (larray[0], 'SR1', 'None'))
            reactants = larray[1].split('|----|')
            products = larray[2].split('|----|')
            catalysts = larray[3].split('|----|')
            solvents = larray[4].split('|----|')
            if larray[1] != '':
                for reactant in reactants:
                    compound = reactant.split('|---|')
                    modelcompounds.append((compound[0], 'SR1'))
                    cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                (larray[0], compound[0], 0, 1, 0))
                    if len(compound) == 2:
                        allcompounds.append((compound[0], compound[1], compartment,  'None'))
                    elif len(compound) < 2:
                        print ('WARNING: Issue with name and ID {} {} reactant'.format(compound, larray[0]))
                        allcompounds.append((compound[0], 'None', compartment,  'None'))
                conn.commit()
            if larray[2] != '':
                for product in products:
                    compound = product.split('|---|')
                    modelcompounds.append((compound[0], 'SR1'))
                    cnx.execute("INSERT INTO reaction_compound VALUES (?,?,?,?,?)",
                                (larray[0], compound[0], 1, 1, 0))
                    if len(compound) == 2:
                        allcompounds.append((compound[0], compound[1], compartment,  'None'))
                    elif len(compound) < 2:
                        print ('WARNING: Issue with name and ID {} {} product'.format(compound, larray[0]))
                        allcompounds.append((compound[0], 'None', compartment, 'None'))
                conn.commit()

            if larray[3] != '':
                for catalyst in catalysts:
                    compound = catalyst.split('|---|')
                    if len(compound) == 2:
                        cnx.execute("INSERT INTO reaction_catalysts VALUES (?,?,?)",
                                    (larray[0], compound[0], compound[1]))
                    elif len(compound) < 2:
                        print ('WARNING: Issue with name and ID {} {} catalyst'.format(compound, larray[0]))
                        cnx.execute("INSERT INTO reaction_catalysts VALUES (?,?,?)",
                                    (larray[0], compound[0], 'None'))
                conn.commit()
            if larray[4] != '':
                for solvent in solvents:
                    compound = solvent.split('|---|')
                    if len(compound) == 2:
                        cnx.execute("INSERT INTO reaction_solvents VALUES (?,?,?)",
                                    (larray[0], compound[0], compound[1]))
                    elif len(compound) < 2:
                        print ('WARNING: Issue with name and ID {} {} solvent'.format(compound, larray[0]))
                        cnx.execute("INSERT INTO reaction_solvents VALUES (?,?,?)",
                                    (larray[0], compound[0], 'None'))
                conn.commit()
            cnx.execute("INSERT INTO reaction_spresi_info VALUES (?,?,?,?,?,?)",
                        (str(larray[0]), str(larray[5]), str(larray[6]), str(larray[7]), str(larray[8]),
                         larray[9].decode('utf-8')))
            conn.commit()

    allcompounds = list(set(allcompounds))
    modelcompounds = list(set(modelcompounds))
    cnx.executemany("INSERT INTO model_compound VALUES (?,?)", modelcompounds)
    cnx.executemany("INSERT INTO compound VALUES (?,?,?,?)", allcompounds)
    conn.commit()
