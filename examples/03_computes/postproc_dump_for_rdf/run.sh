#!/bin/bash
# Compute RDF from the production dump using the external pyrdf utility.
python ../../../utilities/LmpXrd/pyxrd/pyrdf.py in.data_init_mod.dat ../o.dump.prod.lammpstrj 4.0 0.01 10 1 -atomtype='atomic' -rdffile='o.rdf'
