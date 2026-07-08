#!/bin/bash
path_lmpxrd="../../../utilities/LmpXrd/"
path_data="../"

# compute rdf
python $path_lmpxrd/pyxrd/pyrdf.py $path_data/in.data $path_data/o.dump.all.lammpstrj 25.0 0.1 1 1 -atomtype='full' -rdffile='o.rdf_py' > log.rdf

# compute xrd
python $path_lmpxrd/pyxrd/pyxrd.py $path_data/in.data o.rdf_py.dat form_factors.dat 25.0 0.1 20.0 0.01 -atomtype='full' > log.xrd
