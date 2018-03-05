Parameter Documentation
=======================
.. toctree::
   :maxdepth: 2

RetSynth's primary script is: rs.py (main script)

RetSynth Parameters
---------------------

``--targets: List of target molecules and the selected organisms``

``--output_path: Destination for output files``

``--python_glpk_connector_package: python package to use to connect to glpk solver software (can use PYGLPK (GLPK) note that pyglpk only works with glpk package 4.3 or lower or PULP (default))``

``--solver_time_limit: time limit for solver to solve shortest path (note: only function with PULP python solver package)``

``--inchidb: Retrieve InChI's and use them as compound IDs in the database``

``--generate_database: Generates metabolic database``

``--database: Saved metabolic database (.db file)``

``--generate_database_constraints: Generates full stoichometric matrix for a metabolic database`` 

``--database_constraints: Saved stoichmetric matrix for a metabolic database (.constraints file)``

``--kbase: Set to build database with Kbase data``

``--kbase_dump_directory: Directory of metabolic networks used to construct metabolic database (xml files)``
``--metacyc: Set to build database with metacyc data``

``--metacyc_addition: SBML file of reactions in the metacyc database to be integrated in to a metabolic database``

``--metacyc_reaction_type: Type of reactions being added to database from metacyc files bio (default) or chem``
``--translation_file: Translation file to connect metacyc and kbase IDs``

``--SPRESI: Set to build database with SPRESI information``

``--spresi_dump_directory: Directory of SPRESI files (.rdf)``

``--spresi_reaction_type: Type of reactions being added to database from spresi files chem (default) or bio``

``--mine: Set to build database with mine data``

``--mine_dump_directory: Directory of MINE files (.msp)``

``--mine_reaction_type: Type of reactions being added to database from mine files bio (default) or chem``

``--atlas: Set to build database with atlas data``

``--atlas_dump_directory: Directory of ATLAS files (.csv)``

``--atlas_reaction_type: Type of reactions being added to database from atlas files bio (default) or chem``

``--kegg: Set to build database with Kegg data``

``--kegg_organism_type: Define type of organisms reactions to be added to the database bacteria (default) or archea``

``--kegg_number_of_organisms: Define number of organisms from kegg database to add``

``--kegg_number_of_organism_pathways: Define number of pathways from an organism to add``

``--kegg_reaction_type: ype of reactions being added to database from atlas files bio (default) or chem``

``--flux_balance_analysis: Run flux balance analysis``

``--media_for_FBA: Run flux balance analysis on a glucose media``

``--knockouts: Simulates metabolism using FBA for each reaction knockout``

``--limit_reactions: Limit the number of reactions in a identified pathway``

``--limit_cycles: Limits the number of cycle checks``

``--evaluate_reactions: Defines which type of reactions (bio, chem, or all (default)) to be evaluated in identifying shortest paths``

``--start_compounds: List of compounds that the shortest path should start from (use in place of target organism)``

``--k_number_of_paths: Specifies the number of shortest paths wanted by the user``

``--multiple_solutions: Find all shortest paths (True (defualt) or False)``

``--cycles: Elimate cycles when finding shortest paths (True (default) or False)``

``--tanimoto_threshold: Set the tanimoto threshold to identify related compounds``

``--processors: Number of processors to use when solving for shortest path (default 4)``

``--figures: generates figures``

``--images: Set whether to use chemical images or round nodes in output figures (default is True, meaining images will be used)``


