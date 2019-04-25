[![Build Status](https://travis-ci.org/sandialabs/RetSynth.svg?branch=master)](https://travis-ci.org/sandialabs/RetSynth)

# RetSynth

The overaching goal of RetSynth is to streamline the arduous and complex step of selecting enzyme/reactions pairs to produce a target compound for bioengineering microbial organisms. 

## Documentation

See documentation at http://sandialabs.github.io/RetSynth/

## Build
To download RetSynth please follow the instructions below.  The difficult part of ensuring that RetSynth can run is installing the non-python dependencies which include (THESE NEED TO BE INSTALLED FIRST PRIOR TO PYTHON LIBRARIES BELOW):
	
    GNU/GLPK (v4.64) 	 Download from the website http://ftp.gnu.org/gnu/glpk/
	    
    GraphViz (only needed if --figures_graphviz is run)     Download from the website http://graphviz.org/ or using MacPorts, needs Graphviz version 2.38 or lower          

Additionally you will need the python packages glpk, cobra, pygraphviz (only if you want to use Graphviz to generate figures), beautifulsoup and python-libsmbl.

```bash
git clone https://github.com/sandialabs/RetSynth.git
pip install -r requirements.txt
python setup.py install
```

### Accessing GUI
-------------

RetSynth is also distributed as a Windows and Linux/MacOSX GUI - these are separate branches on GitHub.

Accessing Windows Version:
* https://github.com/sandialabs/RetSynth/tree/RetSynthGUIwindows

Accessing MacOSX and Linux Versions:
* https://github.com/sandialabs/RetSynth/tree/RetSynthGUImac_linux


### Dependencies
-------------
RetSynth is currently only tested to work under Python 2.6 and 2.7.

* pulp==1.6.8, only tested with glpk version 4.64
* cobra==0.14.1
* bs4
* pygraphviz==1.3.1
* python-libsbml-experimental==5.10.0

## Publication (please cite)

## Other python dependencies (references) integrated into RetSynth
1. <a href="https://dx.doi.org/doi:10.1093/bioinformatics/btx185" target="_blank">Mackinac: a bridge between ModelSEED and COBRApy to generate and analyze genome-scale metabolic models</a>
2. <a href="http://dx.doi.org/doi:10.1093/nar/gkt1099" target="_blank">PATRIC, the bacterial bioinformatics database and analysis resource</a>
3. <a href="http://dx.doi.org/doi:10.1186/1752-0509-7-74" target="_blank">COBRApy: COnstraints-Based Reconstruction and Analysis for Python</a>

## License

BSD - 3-Clause Copyright 2017 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
