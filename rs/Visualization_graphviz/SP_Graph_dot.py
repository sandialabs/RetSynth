from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Functions to generate figures'

import os
import re
import glob
import httplib
import urllib2
import warnings
from numpy import linspace
from copy import deepcopy
from sys import platform
import numpy as np
import pubchempy
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
from matplotlib import colors
from PIL import Image
import pylab as pl
pl.switch_backend('agg')
from Visualization_graphviz import chempyutilgraph_edit as cuge
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_renderer
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_renderer
elif platform == "win32" or platform == 'win64':
    raise ImportError('Cannot translate RDF file on windows machine')

warnings.filterwarnings("ignore")
edgeweights = []
edgeweights.append([0, .99999999])
edgeweights.append([1, 99.99999999])
edgeweights.append([100, 199.99999999])
edgeweights.append([200, 299.99999999])
edgeweights.append([300, 399.99999999])
edgeweights.append([400, 499.99999999])
edgeweights.append([500, 599.99999999])
edgeweights.append([600, 699.99999999])
edgeweights.append([700, 799.99999999])
edgeweights.append([800, 899.99999999])
edgeweights.append([900, 1000.00000])
PATH = os.path.dirname(os.path.abspath(__file__))
start = 0.0
stop = 0.4
number_of_lines = 10
cm_subsection = linspace(start, stop, number_of_lines)
colors_range = [cm.Reds(x) for x in cm_subsection]
colors_range = [re.sub('\(|\)', '', str(x)) for x in colors_range]
colors_range.append('0.0, 0.0, 0.0, 1.0')
colors_range.reverse()

class GraphDot(object):
    """
    Generates figure of Shortest Path (external reactions and compounds needed to produce
    target compound)
    """
    def __init__(self, db, output_path, incpds, inrxns, temp_img_path, FBA=False):
        self.DB = db
        conn, cnx = self.DB.connect_to_database()
        Q = cnx.execute("""SELECT DISTINCT compartment FROM compound""")
        hits = Q.fetchall()
        self.compartments = [i[0] for i in hits]
        conn.close()
        output_path = re.sub('/$', '', output_path)
        try:
            os.mkdir(output_path+'/solution_figures/')
        except OSError:
            pass
        self.temp_imgs_path = temp_img_path
        self.GRAPHPATH = output_path+'/solution_figures'
        self.output_path = output_path
        self.incpds = incpds
        self.inrxns = inrxns
        self.FBA = FBA
        self.IN = indigo.Indigo()
        self.IR = indigo_renderer.IndigoRenderer(self.IN)
        if self.FBA is not False:
            self.FBA_weights()
            self.generate_flux_legend()

    def generate_flux_legend(self):
        '''Generate legend for flux figures'''
        start_l = 0.5
        stop_l = 0.9
        number_of_lines = 10
        cm_subsection_l = linspace(stop_l, start_l, number_of_lines)
        colors_range_l = [cm.Reds(x) for x in cm_subsection_l]
        colors_range_l.insert(0, (0.0, 0.0, 0.0, 1.0))        
        cmap = colors.ListedColormap(colors_range_l)
        a = np.array([[0, 1.1]])
        pl.figure(figsize=(9, 9))
        img = pl.imshow(a, cmap=cmap)
        pl.gca().set_visible(False)
        cax = pl.axes([0.1, 0.2, 0.8, 0.1])
        cb = pl.colorbar(orientation="horizontal", label='Flux through a reaction',
                         ticks=[.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1], cax=cax)
        cb.ax.set_xticklabels(["0-.99", "5-99.99", "100-199.99", "200-299.99", "300-399.99",
                               "400-499.99", "500-599.99", "600-699.99", "700-799.99",
                               "800-899.99", "900-1000"], rotation=90)
        pl.savefig(self.output_path+"/solution_figures/colorbarforreactionflux.png")

    def get_mol_from_name(self, cpd_name):
        '''Get molecule for a compound from its' name'''
        try:
            pcp_cpds = pubchempy.get_compounds(cpd_name, 'name')
            if pcp_cpds:
                pcp_cpd = pcp_cpds[0]
                SMILE = pcp_cpd.canonical_smiles
                mol = self.IN.loadMolecule(SMILE)
                return mol
        except (pubchempy.PubChemHTTPError, httplib.BadStatusLine,
                urllib2.URLError, ValueError):
            return None

    def crop_figure(self, path_to_image):
        im = Image.open(path_to_image)
        imageBox = im.getbbox()
        cropped = im.crop(imageBox)
        return cropped

    def alter_name_length(self, path_to_figure, cpdname):
        '''Shorten compound name if it is too long'''
        if len(path_to_figure) > 250:
            remove_variable = len(path_to_figure) - 250
            cpdname = cpdname[:-remove_variable]
        return cpdname

    def get_figure(self, cpdID, cpdname, rxn, type_node):
        '''Get smiles for a compounds'''
        if cpdID.startswith('InChI'):
            cpdID = re.sub('_\w{1}\d{1}$', '', cpdID)
            mol = self.IN.loadMolecule(cpdID)
        else:
            for comp in self.compartments:
                if cpdname.endswith('_'+comp):
                    cpdname = re.sub('_'+comp+'$', '', cpdname)
            cpdname = re.sub('_', ' ', cpdname)
            mol = self.get_mol_from_name(cpdname)
            if not mol:
                cpd_name = re.sub(' ', '-', cpdname)
                mol = self.get_mol_from_name(cpd_name)
                if not mol:
                    cpd_name_p = cpd_name+'+'
                    mol = self.get_mol_from_name(cpd_name_p)
                    if not mol:
                        cpd_name_m = cpd_name+'-'
                        mol = self.get_mol_from_name(cpd_name_m)
        if mol:
            if type_node == 'synthetic':
                self.IN.setOption("render-base-color", "0, 0, 1")
            elif type_node == 'internal':
                self.IN.setOption("render-base-color", "0.3, 0.3, 0.3")
            else:
                self.IN.setOption("render-base-color", "1, 0, 0")
            self.IN.setOption("render-relative-thickness", "2")
            self.IN.setOption("render-image-size", "300,200")
            self.IN.setOption("render-margins", "40, 0, 0, 0")
            cpdname = self.reformat_inchi(cpdname)
            cpdname = self.alter_name_length(self.temp_imgs_path+'/compound_'+cpdname+'_'+rxn+'.png', cpdname)
            self.IN.setOption("render-output-format","png")
            self.IR.renderToFile(mol, self.temp_imgs_path+'/compound_'+cpdname+'_'+rxn+'.png')
            cropped = self.crop_figure(self.temp_imgs_path+'/compound_'+cpdname+'_'+rxn+'.png')
            cropped.save(self.temp_imgs_path+'/compound_'+cpdname+'_'+rxn+'_cropped.png', transparent=True)
            return True
        else:
            return False

    def reformat_inchi(self, inchi):
        '''Reformates inchi string'''
        inchi = re.sub('/', '_', inchi)
        inchi = re.sub('-', '_', inchi)
        inchi = re.sub(' ', '_', inchi)
        return inchi

    def synthetic_compound_attr(self, cpdID, name, rxn):
        '''
        Specifies attributes for nodes representing external compounds
        '''
        origname = deepcopy(name)
        if self.images is True:
            name = re.sub('_', '-', name)
            namereformat = self.reformat_inchi(name)
            figure_bool = self.get_figure(cpdID, namereformat, rxn, 'synthetic')
            namereformat = self.alter_name_length(self.temp_imgs_path+'/compound_'+namereformat+'_'+rxn+'.png', namereformat)
            if figure_bool:
                self.outputfile_dot.append('    "{}" [color="{}", image="{}", shape={}, label=""];\n'.format(origname, 'None', self.temp_imgs_path+'/compound_'+namereformat+'_'+rxn+'_cropped.png', 'None'))
            else:
                self.outputfile_dot.append('    "{}" [fillcolor={}, color="{}", height={}, width={}, fontsize={}, fontname="{}", label="{}", fontcolor={}];\n'.format(origname, 'None', 'None', '.3', '.3', '12', 'times', name, 'blue'))
        else:
            self.outputfile_dot.append('    "{}" [fillcolor={}, height={}, width={}, fontsize={}, fontname="{}", xlabel="{}", label=""];\n'.format(origname, 'blue', '.3', '.3', '12', 'times', name))
        self.nodes.add(name)

    def internal_compound_attr(self, cpdID, origname, name, rxn):
        '''
        Specifies attributes for nodes representing internal (native) compounds
        '''
        if self.images is True:
            name = re.sub('_', '-', name)
            namereformat = self.reformat_inchi(name)
            figure_bool = self.get_figure(cpdID, namereformat, rxn, 'internal')
            namereformat = self.alter_name_length(self.temp_imgs_path+'/compound_'+namereformat+'_'+rxn+'.png', namereformat)
            if figure_bool:
                self.outputfile_dot.append('    "{}" [color={}, image="{}", shape={}, label=""];\n'.format(origname, 'None', self.temp_imgs_path+'/compound_'+namereformat+'_'+rxn+'_cropped.png', 'None'))
            else:
                self.outputfile_dot.append('    "{}" [fillcolor={}, color={}, height={}, width={}, fontsize={}, fontname="{}", label="{}", fontcolor={}];\n'.format(origname, 'None', 'None', '.2', '.2', '10', 'times', name, 'gray33'))
        else:
            self.outputfile_dot.append('    "{}" [fillcolor={}, height={}, width={}, fontsize={}, fontname="{}", xlabel="{}", label=""];\n'.format(origname, 'gray', '.2', '.2', '10', 'times', name))
        self.nodes.add(origname)

    def add_reactants(self, reactants, rxn):
        '''
        Add nodes to graph for the reactants of reactions
        '''
        for react in reactants:
            reactname = reactants[react]
            if (reactname == 'None'
                    or reactname == 'none' or reactname.startswith('None')):
                reactname = react
            reactname = re.sub('-', '_', reactname)
            reactname = re.sub('_\w{1}0$', '', reactname)
            if react not in self.incpds:
                if self.FBA is not False and self.flux_values[rxn] > 0:
                    countkey = max(self.rxncolors[rxn].keys())
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", penwidth={}, label=""];\n'.format(reactname, rxn, self.rxncolors[rxn][countkey], 3))
                elif self.FBA is not False and self.flux_values[rxn] == 0:
                    countkey = max(self.rxncolors[rxn].keys())
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", penwidth={}, label=""];\n'.format(reactname, rxn, self.rxncolors[rxn][countkey], 3))
                else:
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", label=""];\n'.format(reactname, rxn, 'black'))
                if reactname not in self.nodes:
                    self.synthetic_compound_attr(react, reactname, rxn)

            else:
                self.count += 1
                if self.FBA is not False:
                    countkey = max(self.rxncolors[rxn].keys())
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", penwidth={}, label=""];\n'.format(reactname+'_'+str(self.count), rxn, self.rxncolors[rxn][countkey], 3))
                else:
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", label=""];\n'.format(reactname+'_'+str(self.count), rxn, 'black'))
                self.internal_compound_attr(react, reactname+'_'+str(self.count), reactname, rxn)

    def add_products(self, products, rxn):
        '''
        Add nodes to graph for the products of reactions
        '''
        for prod in products:
            prodname = products[prod]
            if prodname == 'None' or prodname == 'none' or prodname.startswith('None'):
                prodname = prod
            prodname = re.sub('-', '_', prodname)
            prodname = re.sub('_\w{1}0$', '', prodname)
            prodnamereformat = self.reformat_inchi(prodname)
            if prod not in self.incpds:
                if self.FBA is not False and self.flux_values[rxn] > 0:
                    countkey = max(self.rxncolors[rxn].keys())
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", penwidth={}, label=""];\n'.format(rxn, prodname, self.rxncolors[rxn][countkey], 3))
                elif self.FBA is not False and self.flux_values[rxn] == 0:
                    countkey = max(self.rxncolors[rxn].keys())
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", penwidth={}, label=""];\n'.format(rxn, prodname, self.rxncolors[rxn][countkey], 3))
                else:
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", label=""];\n'.format(rxn, prodname, 'black'))
                if prod == self.target:
                    if self.images is True:
                        figure_bool = self.get_figure(prod, prodname, rxn, 'target')
                        prodnamereformat = self.alter_name_length(self.temp_imgs_path+'/compound_'+prodnamereformat+'_'+rxn+'.png', prodnamereformat)
                        if figure_bool:
                            self.outputfile_dot.append('    "{}" [color={}, image="{}", shape={}, label=""];\n'.format(prodname, 'None', self.temp_imgs_path+'/compound_'+prodnamereformat+'_'+rxn+'_cropped.png', 'None'))
                        else:
                            self.outputfile_dot.append('    "{}" [fillcolor={}, color="{}", height={}, width={}, fontsize={}, fontname="{}", label="{}", fontcolor={}];\n'.format(prodname, 'None', 'None', '.4', '.4', '12', 'times', prodname, 'red'))
                    else:
                        self.outputfile_dot.append('    "{}" [fillcolor={}, height={}, width={}, fontsize={}, fontname="{}", xlabel="{}", label=""];\n'.format(prodname, 'red', '.4', '.4', '12', 'times', prodname))
                    self.nodes.add(prodname)
                elif prodname not in self.nodes:
                    self.synthetic_compound_attr(prod, prodname, rxn)
            else:
                self.count += 1
                if self.FBA is not False:
                    countkey = max(self.rxncolors[rxn].keys())
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", penwidth={}, label=""];\n'.format(rxn, prodname+'_'+str(self.count), self.rxncolors[rxn][countkey], 3))
                else:
                    self.outputfile_dot.append('    "{}" -> "{}"    [color="{}", label=""];\n'.format(rxn, prodname+'_'+str(self.count), 'black'))
                self.internal_compound_attr(prod, prodname+'_'+str(self.count), prodname, rxn)

    def FBA_weights(self):
        '''
        If FBA was run alter the thickness of the edges between
        reactions and compounds to represent the reaction activityhttp://localhost:8888/notebooks/software/retrosynth/rs/rsjuypter.ipynb#
        '''
        self.rxncolors = {}
        self.flux_values = {}
        self.temp = []
        for r in self.FBA.index:
            for count, e in enumerate(edgeweights):
                if (abs(float(round(self.FBA[r], 2))) >= float(e[0]) and
                        abs(float(round(self.FBA[r], 2))) <= float(e[1])):
                    self.rxncolors[r] = {}
                    self.rxncolors[r][count] = colors_range[count]
                    self.flux_values[r] = round(abs(float(self.FBA[r])), 2)

    def sc_graph(self, target_compound_ID, target_organism_ID, paths, images=True, RP=None):
        '''
        Build and output graph to a file, reaction nodes added in this function
        '''
        self.outputfile_dot = []
        self.target = target_compound_ID
        self.images = images
        self.outputfile_dot.append('digraph None{\n')
        self.outputfile_dot.append('    graph [bb="0,0,186,520", landscape=false, margin="2.5,.7", ranksep=.7];\n')
        self.outputfile_dot.append('    node [fixedsize=true, fontcolor=black, label="\N", shape=circle, style=filled];\n')
        self.count = 0
        self.store_rxns = set()
        self.nodes = set()

        def load_into_dot_file(rxns):
            for rxn in rxns:
                name = rxns[rxn]['name']
                if rxn not in self.store_rxns:
                    org = rxns[rxn]['organisms'][0]
                    name = re.sub('_\w{1}0$', '', name)
                    if (self.DB.get_proteins(rxn, rxns[rxn]['organisms'][0]) != 'None') and (name == 'None' or name.startswith('RXN') or name.startswith('rxn') or name.endswith('RXN')):
                        xlabel = self.DB.get_proteins(rxn, org)
                    else:
                        if self.DB.get_proteins(rxn, rxns[rxn]['organisms'][0]) != 'None':
                            xlabel = name+'\t'+self.DB.get_proteins(rxn, org)
                        elif self.DB.get_genes(rxn, rxns[rxn]['organisms'][0]) != 'None':
                            xlabel = self.DB.get_genes(rxn, rxns[rxn]['organisms'][0])
                        else:
                            xlabel = rxn
                    if rxn not in self.inrxns:
                        self.outputfile_dot.append('    "{}" [fillcolor={}, fontname="{}", fontsize={}, shape={}, height={}, width={}, xlabel="{}", label=""];\n'.format(rxn, 'black', 'times', '12', 'rect', '.3', '.15', xlabel))
                    else:
                        self.outputfile_dot.append('    "{}" [fillcolor={}, fontname="{}", fontsize={}, shape={}, height={}, width={}, xlabel="{}", label=""];\n'.format(rxn, 'gray', 'times', '12', 'rect', '.3', '.15', xlabel))                        
                    self.store_rxns.add(rxn)
                    if rxns[rxn]['direction'] == 'forward':
                        self.add_reactants(rxns[rxn]['reactants'], rxn)
                        self.add_products(rxns[rxn]['products'], rxn)
                    elif rxns[rxn]['direction'] == 'reverse':
                        self.add_reactants(rxns[rxn]['products'], rxn)
                        self.add_products(rxns[rxn]['reactants'], rxn)
        if RP:
            for path_count, rxns in paths.iteritems():
                if path_count in RP.correctly_typed_paths:
                    load_into_dot_file(rxns)
        else:
           for path_count, rxns in paths.iteritems():
               load_into_dot_file(rxns)

        self.outputfile_dot.append('    subgraph None{\n\tlabel="test"\n\t}\n')
        self.outputfile_dot.append('}\n')
        cpdname = self.DB.get_compound_name(target_compound_ID)
        if cpdname.startswith('None'):
            cpdname = self.reformat_inchi(target_compound_ID)
        else:
            cpdname = self.reformat_inchi(cpdname)
        cuge.dot2graph(self.outputfile_dot, 'SC_graph_'+str(cpdname)+'_'+str(self.DB.get_organism_name(str(target_organism_ID)))+'.png', self.GRAPHPATH)
        graph_files = glob.glob(os.path.join(self.GRAPHPATH, '*'))
