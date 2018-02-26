[![Build Status](https://travis-ci.org/sandialabs/RetroSynth.svg?branch=master)](https://travis-ci.org/sandialabs/RetSynth)

# RetSynth

The overaching goal of RetSynth is to streamline the arduous and complex step of selecting enzyme/reactions pairs to produce a target compound for bioengineering microbial organisms. 

## Documentation

See documentation at https://sandialabs.github.io/RetSynth/#

## Build

The difficult part of ensuring that RetSynth can run is installing the non-python dependencies which include:
	
    GNU/GLPK 	 Download from the website http://ftp.gnu.org/gnu/glpk/
    
    GraphViz     Download from the website http://graphviz.org/ or using MacPorts

Additionally you will need the python packages glpk, cobra, pygraphviz, beautifulsoup and python-libsmbl.

```bash
git clone https://github.com/sandialabs/RetSynth.git
pip install -r requirements.txt
python setup.py install
```

### Dependencies
-------------
RetSynth is currently only tested to work under Python 2.6 and 2.7.

* glpk==0.3 or pulp==1.6.8 (preferred)
* cobra==0.6.2
* beautifulsoup==4-4.4.1-3
* pygraphviz==1.3.1
* python-libsbml-experimental==5.10.0

## Publication (please cite)

### License
-----------
BSD - 3-Clause Copyright 2017 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
=======
