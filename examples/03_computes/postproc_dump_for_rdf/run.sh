#!/bin/bash
# Compute RDF from the production dump using the external pyrdf utility.
python ../../../utilities/LmpXrd/pyxrd/pyrdf.py in.data_init_mod.dat ../o.dump.prod.lammpstrj 4.0 0.01 1 10 -atomtype='atomic' -rdffile='o.rdf'
