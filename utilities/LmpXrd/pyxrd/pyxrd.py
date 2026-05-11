###############################################################################
# MIT License
#
# Copyright (c) 2024 ArisSgouros
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

import os
import sys
import ast
import numpy as np
import math as math
import argparse

kTol = 1e-6
kRadToDegree = 180.0/np.pi

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('datafile', type=str, help='Path of the lammps data file')
parser.add_argument('rdffile', type=str, default="", help='Path of the RDF file')
parser.add_argument('formfactorfile', type=str, default="FORM_FACTORS.dat", help='Path of the form factors')
parser.add_argument('rmax', type=float, help='Max radious')
parser.add_argument('dr', type=float, help='Step interval')
parser.add_argument('qmax', type=float, help='Max wavevector')
parser.add_argument('dq', type=float, help='wavevector interval')
parser.add_argument('-atomtype', type=str, default='full', help='Set the atomtype in the Lammps data file (e.g., full, atomic, etc.)')
parser.add_argument('-lambda_inc', type=float, default=1.5406, help='wavelength of incident beam [default: 1.5406 AA, Cu Ka ]')

def FormFact(p, qq):
   fa = p["c"]
   for ii in range(4):
      fa += p["a"][ii] * np.exp(-p["b"][ii]*pow(qq/(4*np.pi),2))
   return fa

if __name__ == "__main__":
   args = parser.parse_args()
   file_data = args.datafile
   file_rdf  = args.rdffile
   file_FF   = args.formfactorfile
   rmax         = args.rmax
   dr           = args.dr
   qmax         = args.qmax
   dq           = args.dq
   atomtype     = args.atomtype
   lambda_inc   = args.lambda_inc

   print( "*Parameters of the calculation*")
   print()
   print( "file_data   :",file_data)
   print( "file_rdf    :",file_rdf)
   print( "rmax           :",rmax)
   print( "dr             :",dr)
   print( "qmax           :",qmax)
   print( "dq             :",dq)

   # Deal with the format of the Lammps data file
   if atomtype == 'full':
      col_type = 2
   elif atomtype == 'atomic':
      col_type = 1
   else:
      print('unsupported format of Lammps data files for the atom type ' + atomtype)
      sys.exit()

   # Initialize the binning
   rbins = int(rmax / dr)
   r_range = [(ii+0.5)*dr for ii in range(rbins)]
   qbins = int(qmax / dq)
   q_range = [(ii+0.5)*dq for ii in range(qbins)]

   # Convert qq to 2*theta
   two_theta_range = []
   for qq in q_range:
      two_theta = -1.0
      sintheta = lambda_inc*qq/4.0/np.pi
      if abs(sintheta) < 1.0:
         two_theta = 2.0*np.arcsin(sintheta)*kRadToDegree
      two_theta_range.append(two_theta)

   print( "Parsing the atom types from the Lammps datafile:",file_data,"..")

   lines = []
   iline = 0

   with open(file_data,"r") as openFileObject:
      for cur_line in openFileObject:

         if "atom types" in cur_line:
            n_atom_type = int(cur_line.split()[0])
         if "atoms" in cur_line:
            n_atom = int(cur_line.split()[0])
         if "Masses" in cur_line:
            line_start_type = iline + 2
            line_end_type = line_start_type + n_atom_type
         if "Atoms" in cur_line:
            line_start_atom = iline + 2
            line_end_atom = line_start_atom + n_atom

         lines.append(cur_line)
         iline += 1

   type_of_species = {}
   species_of_type = {}
   for cur_line in range(line_start_type,line_end_type):
      cur_lineSplit = lines[cur_line].split()
      species = cur_lineSplit[3]
      type = int(cur_lineSplit[0])

      if species in type_of_species.keys():
         type_of_species[species].append(type)
      else:
         type_of_species[species] = [type]

      species_of_type[type] = species

   for species in type_of_species:
      print( species, type_of_species[species])

   ids_of_species = {}
   for species in type_of_species:
      ids_of_species[species] = []

   for cur_line in range(line_start_atom,line_end_atom):
      cur_lineSplit = lines[cur_line].split()
      id = int(cur_lineSplit[0])
      type = int(cur_lineSplit[col_type])
      species = species_of_type[type]
      ids_of_species[species].append(id)

   N = {}
   C = {}
   N["all"] = 0
   C["all"] = 1.0
   for species in type_of_species:
      N[species] = len(ids_of_species[species])
      N["all"] += N[species]

   for species in type_of_species:
      C[species] = float(N[species]) / float(N["all"])
      print( "There are",N[species]," ",species,"atoms (",C[species],"%)")

   print( "There are",N["all"], "atoms in total..")

   #Write the atom IDs for each species in a file
   print( "Printing the atom IDs for each species..")
   for species in type_of_species:
      f = open('o.IDS.'+species+'.dat', 'w')
      f.write(species+"\n")
      for id in ids_of_species[species]:
         f.write(str(id)+" ")
   f.close()

   print( "Reading the partial gr from file",file_rdf,"..")
   gab_all = {}
   f = open(file_rdf,"r")
   f.readline()
   rho = float(f.readline().split()[2])
   print( rho)
   types = f.readline().split()
   types.pop(0)
   for type in types:
      gab_all[type] = [0] * rbins
   f.readline() # CA
   f.readline() # CB
   for ir in range(rbins):
      current_line = f.readline().split()
      current_line.pop(0)
      for it in range(len(types)):
         type = types[it]
         gab_all[type][ir] = float(current_line[it])

   #
   # This section deals with the atomic form factors
   #
   # The structure factors were retrieved from :
   # http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
   print( "Reading the atomic form factors from",file_FF,"..\n")

   FFp = {}
   for species in type_of_species:
      FFp[species] = {"a" : [0.0]*4, "b" : [0.0]*4, "c" : 0.0}

   with open(file_FF,"r") as openFileObject:
      for cur_line in openFileObject:
         cur_line = cur_line.split()
         for species in type_of_species:
            if species == cur_line[0]:
               FFp[species]["a"][0] = float(cur_line[1])
               FFp[species]["a"][1] = float(cur_line[3])
               FFp[species]["a"][2] = float(cur_line[5])
               FFp[species]["a"][3] = float(cur_line[7])
               FFp[species]["b"][0] = float(cur_line[2])
               FFp[species]["b"][1] = float(cur_line[4])
               FFp[species]["b"][2] = float(cur_line[6])
               FFp[species]["b"][3] = float(cur_line[8])
               FFp[species]["c"]    = float(cur_line[9])

   print( "Atomic Form Factors")

   for species in FFp:
      print( "species :", species)
      print( "a", FFp[species]["a"])
      print( "b", FFp[species]["b"])
      print( "c", FFp[species]["c"])
   print()

   print( "Computation of the Faber-Ziman structure factor..\n")

   SFZ_all = {}
   for species_A in type_of_species:
      for species_B in type_of_species:

         SFZ = [0.0] * qbins
         #print( "SFZ computation of",species_A,"->",species_B,"..")

         type = species_A+"_"+species_B
         gab = gab_all[type]
         for iq in range(qbins):
            qq = q_range[iq]
            r_integral = 0.0
            for ir in range(rbins):
               rr = r_range[ir]
               r_integral += dr * pow(rr,2) * math.sin(qq*rr)/(qq*rr) * (gab[ir]-1)

            SFZ[iq] = 1.0 + 4 * np.pi * rho * r_integral
         SFZ_all[type] = SFZ

   g = open('o.All_SFZ.dat', 'w')
   g.write("%-20s" % ("type"))
   for key in gab_all:
      g.write("%-20s" % (key))
   g.write("\n")

   g.write("%-20s" % ("ca"))
   for key in gab_all:
      species_A = key.split("_")[0]
      g.write("%-20s" % (str(C[species_A])))
   g.write("\n")

   g.write("%-20s" % ("cb"))
   for key in gab_all:
      species_B = key.split("_")[1]
      g.write("%-20s" % (str(C[species_B])))
   g.write("\n")

   for iq in range(qbins):
      g.write("%-20.6f" %(q_range[iq]))
      for key in SFZ_all:
         g.write("%-20.6f" % (SFZ_all[key][iq]))
      g.write("\n")
   g.close()

   print( "Computation of the X-ray weighted structure factor..\n")

   Fx = [0.0] * qbins
   for iq in range(qbins):
      qq = q_range[iq]

      denominator = 0.0
      for species_A in type_of_species:
         FA = FormFact(FFp[species_A],qq)
         denominator += C[species_A] * FA
      denominator = pow(denominator,2)

      numerator = 0.0
      for species_A in type_of_species:
         for species_B in type_of_species:
            type = species_A+"_"+species_B
            FA = FormFact(FFp[species_A],qq)
            FB = FormFact(FFp[species_B],qq)

            numerator += FA * FB * C[species_A] * C[species_B] * (SFZ_all[type][iq] - 1.0)

      Fx[iq] = numerator / denominator


   print( "Computation of the coherent scattering intensity..\n")

   Icoh = [0.0] * qbins
   for iq in range(qbins):
      qq = q_range[iq]

      denominator = 0.0
      for species_A in type_of_species:
         FA = FormFact(FFp[species_A],qq)
         denominator += C[species_A] * FA
      denominator = pow(denominator,2)

      sum_c_fsq = 0.0
      for species_A in type_of_species:
         FA = FormFact(FFp[species_A],qq)
         sum_c_fsq += C[species_A] * FA * FA

      Icoh[iq] = Fx[iq] * denominator + sum_c_fsq

   g = open('o.Fx_Icoh.dat', 'w')
   g.write("%-19s %-19s %-19s %-19s %-19s\n" % ("bin", "q", "two_theta", "Fx", "Icoh"))
   for iq in range(qbins):
      g.write( "%-19d %-19.6f %-19.6f %-19.6f %-19.6f\n"  %(iq, q_range[iq], two_theta_range[iq], Fx[iq], Icoh[iq]))
   g.close()
   print( "Done!")

   #
   # Direct calculation of the coherent scattering intensity
   #
   S_total = [0.0] * qbins

   for iq in range(qbins):
      qq = q_range[iq]

      term_B = 0.0

      for species_A in type_of_species:
         for species_B in type_of_species:

            r_integral = 0.0
            gab = gab_all[species_A+"_"+species_B]
            for ir in range(rbins):
               rr = r_range[ir]
               r_integral += dr * pow(rr,2) * math.sin(qq*rr)/(qq*rr) * (gab[ir]-1)

            FA = FormFact(FFp[species_A],qq)
            FB = FormFact(FFp[species_B],qq)

            sum = 4 * np.pi * rho * FA * FB * C[species_A] * C[species_B] * r_integral

            term_B += sum

      term_A = 0.0
      for species_A in type_of_species:
         FA = FormFact(FFp[species_A],qq)
         term_A += C[species_A] * FA * FA

      S_total[iq] = term_B + term_A

   g = open('o.S_XRD.dat', 'w')
   g.write("%-19s %-19s %-19s\n" % ("bin", "q", "S"))
   for iq in range(qbins):
      g.write( "%-19d %-19.6f %-19.6f\n"  %(iq, q_range[iq], S_total[iq]))
   g.close()
   print( "Done!")
