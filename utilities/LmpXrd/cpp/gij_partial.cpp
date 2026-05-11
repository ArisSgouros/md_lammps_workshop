/*
 * MIT License
 *
 * Copyright (c) 2024 ArisSgouros
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>

using namespace std;

std::vector<std::string> tokenize(std::string input_string) {
   std::string buf;
   std::stringstream ss(input_string);
   std::vector<std::string> tokens;
   while (ss >> buf)
      tokens.push_back(buf);
   return tokens;
}

struct Atom{
  int id;
  int mol;
  int type;
  float rx;
  float ry;
  float rz;
};

int main(int argc, char** argv){

   // Parse the command line arguments
   if (argc < 8){
      cout<<endl;;
      cout<<"The required formats are the following:"<<endl;
      cout<<endl;
      cout<<"./gij_partial \"r_max\" \"dr\" \"n_every\" \"interactions\" \"datafile\" \"dumpfile\" \"pairfile\" \"formatfile\"  \"rdffile\""<<endl;
      cout<<endl;
      cout<<"with INTERACTIONS being set to either \"INTER\", \"INTRA\" or \"ALL\""<<endl;
      cout<<endl;
      cout<<"example: ./gij_partial 15 0.5 1 INTER CG.dat CG.lammpstrj gij_pairs.dat atomstyle o.rdf"<<endl;
      cout<<endl;
      cout<<"*The format of the Masses section of the data file should be like the following:"<<endl;
      cout<<endl;
      cout<<" Masses"<<endl;
      cout<<endl;
      cout<<" 1    15.9994   # O"<<endl;
      cout<<" 3    12.0107   # C"<<endl;
      cout<<" 2    12.0107   # C"<<endl;
      cout<<" .."<<endl;
      cout<<endl;
      cout<<"*The format of the pairfile should be like the following (eg. for 3 pairs):"<<endl;
      cout<<endl;
      cout<<" 3            number of pairs"<<endl;
      cout<<endl;
      cout<<" C C"<<endl;
      cout<<" C O"<<endl;
      cout<<" O O"<<endl;
      cout<<" .."<<endl;
      cout<<endl;
      cout<<"exiting.."<<endl;
      return 0;
   }

   float r_max            = atof(argv[1]);
   float r_bin            = atof(argv[2]);
   int n_every            = atoi(argv[3]);
   std::string inter_type = argv[4];
   std::string datafile   = argv[5];
   std::string dumpfile   = argv[6];
   std::string pairfile   = argv[7];
   std::string atom_style    = argv[8];
   std::string rdffile    = argv[9];
   bool export_verbose      = true;
   if (argc > 10) export_verbose = (bool)atoi(argv[10]);

   // Deal with file formatting
   int data_col_id   = -1;
   int data_col_mol  = -1;
   int data_col_type = -1;
   int dump_col_id   = -1;
   int dump_col_rx   = -1;
   int dump_col_ry   = -1;
   int dump_col_rz   = -1;

   // set the format of the data file based on the Lammps convention
   if (atom_style == "full") {
      data_col_id = 0;
      data_col_mol = 1;
      data_col_type = 2;
   } else if (atom_style == "atomic") {
      data_col_id = 0;
      data_col_type = 1;
   } else if (atom_style == "angle") {
      data_col_id = 0;
      data_col_mol = 1;
      data_col_type = 2;
   } else {
      std::cout<<"Error: unsupported atom format "<<atom_style<<std::endl;
      return 1;
   }

   // parse the format from the header of the dump file
   {
      ifstream tmp_dfile(dumpfile.c_str(), ifstream::in);
      std::string current_line;
      for (int ii=0; ii<9; ii++) {
         getline(tmp_dfile, current_line); // header of dump file section
      }
      std::vector <std::string> tokens = tokenize(current_line);
      tokens.erase(tokens.begin()); // rmv first
      tokens.erase(tokens.begin()); // rmv second
      for (int icol=0; icol<tokens.size(); icol++) {
         std::string quantity = tokens[icol];
         if (quantity == "id") {
            dump_col_id = icol;
         } else if (quantity == "xu" || quantity == "x" || quantity == "xs") {
            dump_col_rx = icol;
         } else if (quantity == "yu" || quantity == "y" || quantity == "ys") {
            dump_col_ry = icol;
          } else if (quantity == "zu" || quantity == "z" || quantity == "zs") {
            dump_col_rz = icol;
         }
      }
      tmp_dfile.close();
   }

   cout<<endl;
   cout<<"Parameters of the computation:"<<endl;
   cout<<"*r_max: "<<r_max<<endl;
   cout<<"*r_bin: "<<r_bin<<endl;
   cout<<"*the lammps dump file is processed every "<<n_every<<" frames"<<endl;

   int interact_flag = 0;
   if (inter_type=="INTER"){
      cout << "*gij based in INTERMOLECLAR interactions\n";
      interact_flag = 1;
   }else if (inter_type=="INTRA"){
      cout << "*gij based in INTRAMOLECLAR interactions\n";
      interact_flag = 2;
   }else{
      cout << "*gij based in ALL interactions\n";
      interact_flag = 0;
   }

   float r_max2 = r_max * r_max;
   int n_bin    = (int)(r_max / r_bin);
   float ir_bin = 1.0 / r_bin;

   cout<<"\nReading the lammps datafile with name "<<datafile<<"..\n"<<endl;

   int n_atoms=0, n_types=0, iline = 0;
   int atoms_line_start = 0, masses_line_start = 0, atom_types_line = 0, pairs_line_start = 0;

   std::vector <std::string> lines_of_file;
   std::vector <std::string> tokens;
   lines_of_file.clear();
   ifstream data_file(datafile.c_str(), ifstream::in);
   while (data_file.good()) {
      std::string current_line;
      getline(data_file, current_line);
      lines_of_file.push_back(current_line);
      if (current_line.find("atoms") != std::string::npos){
         n_atoms = atoi(current_line.c_str());
         cout<<"number of atoms: "<<n_atoms<<endl;
      }
      if ((current_line.find("atom") != std::string::npos) && (current_line.find("types") != std::string::npos)){
         n_types = atoi(current_line.c_str());
         cout<<"number of types: "<<n_types<<endl;
      }
      if (current_line.find("Atoms") != std::string::npos)
         atoms_line_start = iline+2;
      if (current_line.find("Masses") != std::string::npos)
         masses_line_start = iline+2;
       iline++;
   }
   data_file.close();

   // A dictionary that maps (int) Lammps atom types to (str) species labels
   std::map<int, std::string > species_of_type;

   // Read the masses section and assign species to types
   for (int iline = masses_line_start; iline < masses_line_start + n_types; iline++){
      tokens = tokenize(lines_of_file[iline]);
      int type = stoi(tokens[0])-1;
      std::string species = tokens[3];
      species_of_type[type] = species;
   }

   cout<<"\nSpecies assigned to each lammps atom type:"<<endl;
   for (std::map<int, std::string >::iterator it = species_of_type.begin(); it != species_of_type.end(); it++)
      cout<<it->first<<" "<<it->second<<endl;
   cout<<endl;

   // Create containers that assign number to species and species to numbers
   // This is done to map multiple lammps types to unique species.
   std::map<std::string, int> num_of_species;
   std::vector<std::string> species_of_num;
   int snum = 0;
   for (std::map<int, std::string >::iterator it = species_of_type.begin(); it != species_of_type.end(); it++){
      int type = it->first;
      std::string species;
      species.assign(species_of_type[type]);
      if (num_of_species.count(species) == 0) {
         num_of_species[species] = snum;
         species_of_num.push_back(species);
         cout<<"..assigning num "<<snum<<" to species "<<species_of_num[snum]<<endl;
         snum += 1;
      }
   }
   cout<<endl;

   // Number of unique species
   int n_species = num_of_species.size();

   // A dictionaty that contains the ids for each species
   std::vector<std::vector<int>> ids_of_species;

   std::vector<int> ids;
   for (int ii = 0; ii < n_species; ii++)
      ids_of_species.push_back(ids);

   // A vector array that contains all atom classes
   std::vector<Atom> atoms;
   atoms.resize(n_atoms);

   // Read the atom section
   //std::vector <std::string> lines_of_file;
   //std::vector <std::string> tokens;
   for (int iline = atoms_line_start; iline < atoms_line_start + n_atoms; iline++){
      tokens = tokenize(lines_of_file[iline]);
      Atom atom;
      atom.id   = stoi(tokens[data_col_id])-1;;
      if (data_col_mol >= 0) {
         atom.mol  = stoi(tokens[data_col_mol])-1;
      } else {
         atom.mol = -1;
      }
      atom.type = stoi(tokens[data_col_type])-1;
      atoms[atom.id] = atom;
      int snum = num_of_species[species_of_type[atom.type]];
      ids_of_species[snum].push_back(atom.id);
   }

   // Write the Ids of each species to files
   if (export_verbose) {
      for (int ii=0; ii < n_species; ii++){
         string id_file_name;
         id_file_name += "o.IDS."+species_of_num[ii]+".dat";
         ofstream idfile;
         idfile.open (id_file_name.c_str());
         cout<<"..writting IDs of species "<<species_of_num[ii]<<" to file " <<id_file_name<<endl;
         for (std::vector<int>::iterator id = ids_of_species[ii].begin(); id != ids_of_species[ii].end(); id++)
            idfile<<(*id)<<" ";
         idfile.close();
      }
   }

   // Read the pairs file
   cout<<"\nReading the pairs of species from file: "<<pairfile<<endl;
   ifstream pair_file(pairfile.c_str(), ifstream::in);
   // Define an array of strings to hold the contents of the file.

   lines_of_file.clear();
   iline = 0;
   int n_pairs=0;
   while (pair_file.good()) {
      // A temporary string for the current line of the file.
      std::string current_line;
      getline(pair_file, current_line);
      // Add current line to file's array of lines.
      lines_of_file.push_back(current_line);

      if (current_line.find("pairs") != std::string::npos){
         n_pairs = atoi(current_line.c_str());
         cout<<"number of pairs: "<<n_pairs<<endl;
      }
      if (current_line.find("Pairs") != std::string::npos)
         pairs_line_start = iline+2;
      iline += 1;
   }
   pair_file.close();

   std::vector< std::pair <int,int> > pairs;
   for (int iline = pairs_line_start; iline < pairs_line_start + n_pairs; iline++){
      tokens = tokenize(lines_of_file[iline]);
      if (lines_of_file[iline].empty()){
         cout<<"#error! Problem with reading pairs"<<endl;
         return 1;
      }
      std::string species_A, species_B;
      species_A.assign(tokens[0]);
      species_B.assign(tokens[1]);
      pairs.push_back( make_pair(num_of_species[species_A],num_of_species[species_B]) );
   }

   cout<<"\nlist of pairs to be examined"<<endl;
   for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();
        it_pair != pairs.end(); it_pair++)
      cout<<it_pair->first<<" "<<it_pair->second<<" ( "<<species_of_num[it_pair->first]<<" "<<species_of_num[it_pair->second]<<" )"<<endl;

   // find the number and concentration of each species
   std::vector<float> c_of_species;
   std::vector<int> N_of_species;
   for (int ii=0; ii < n_species; ii++){
      N_of_species.push_back( ids_of_species[ii].size() );
      c_of_species.push_back( (float)N_of_species[ii] / ((float)n_atoms) );
   }

   // find the number of interactions
   cout<<"\nfind the number of interactions.."<<endl;
   float N_inter[n_species][n_species];
   for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();
      it_pair != pairs.end(); it_pair++){
      int species_A = it_pair->first;
      int species_B = it_pair->second;
      N_inter[species_A][species_B] = 0;
      for (std::vector<int>::iterator itA = ids_of_species[species_A].begin(); itA != ids_of_species[species_A].end(); itA++){
         int iatom = (*itA);
         for (std::vector<int>::iterator itB = ids_of_species[species_B].begin(); itB != ids_of_species[species_B].end(); itB++){
            int jatom = (*itB);
            if (atoms[iatom].mol != atoms[jatom].mol){;
               N_inter[species_A][species_B] += 1.0;
            }
         }
      }
      N_inter[species_A][species_B] /= (float)N_of_species[species_A];
      cout<<species_A<<" "<<species_B<<" = "<<N_inter[species_A][species_B]<<endl;
   }

   // Initialize the gij arrays
   int partial_hist_ij[n_species][n_species][n_bin];
   float partial_gij[n_species][n_species][n_bin];
   for (int i = 0; i < n_species; ++i)
      for (int j = 0; j < n_species; ++j)
         for (int k = 0; k < n_bin; ++k){
            partial_hist_ij[i][j][k] = 0;
            partial_gij[i][j][k] = 0.0;
         }

   cout<<"\nComputation of gij from lammps dumpfile: "<<dumpfile<<endl;
   int iframe = 0;
   int n_samples = 0;
   float volume = 0.0;

   {
      ifstream dump_file(dumpfile.c_str(), ifstream::in);

      for (int iframe = 0; ;iframe++){
         // A temporary string for the current line of the file.
         std::string current_line;

         // Check if end of file
         getline(dump_file, current_line); //ITEM: TIMESTEP
         if (not dump_file.good())
            break;

         // Check if this is a frame we need to process.
         if (iframe % n_every != 0){
            for (int ii=0; ii<n_atoms+8; ii++){
               getline(dump_file,current_line);
            }
            continue;
         }
         else
            n_samples += 1;

         float lx, ly, lz, xlo, ylo, zlo, xhi, yhi, zhi;

         getline(dump_file, current_line);
         getline(dump_file, current_line); // ITEM: NUMBER OF ATOMS
         getline(dump_file, current_line);
         getline(dump_file, current_line); // ITEM: BOX BOUNDS pp pp pp
         getline(dump_file, current_line);
         tokens = tokenize(current_line);
         xlo = stof(tokens[0]);
         xhi = stof(tokens[1]);
         getline(dump_file, current_line);
         tokens = tokenize(current_line);
         ylo = stof(tokens[0]);
         yhi = stof(tokens[1]);
         getline(dump_file, current_line);
         tokens = tokenize(current_line);
         zlo = stof(tokens[0]);
         zhi = stof(tokens[1]);

         lx = xhi - xlo;
         ly = yhi - ylo;
         lz = zhi - zlo;
         float ilx = 1.0 / lx;
         float ily = 1.0 / ly;
         float ilz = 1.0 / lz;

         volume += lx * ly * lz;

         getline(dump_file, current_line); // ITEM: ATOMS id mol type xu yu zu

         for (int ii = 0; ii < n_atoms; ii++){
            getline(dump_file, current_line);
            tokens = tokenize(current_line);
            int id = stoi(tokens[dump_col_id]) -1;
            atoms[id].rx   = stof(tokens[dump_col_rx]);
            atoms[id].ry   = stof(tokens[dump_col_ry]);
            atoms[id].rz   = stof(tokens[dump_col_rz]);
         }
         for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();
              it_pair != pairs.end(); it_pair++){
            int species_A = it_pair->first;
            int species_B = it_pair->second;
            for (std::vector<int>::iterator itA = ids_of_species[species_A].begin(); itA != ids_of_species[species_A].end(); itA++){
               int iatom = (*itA);
               float irx = atoms[iatom].rx;
               float iry = atoms[iatom].ry;
               float irz = atoms[iatom].rz;
               for (std::vector<int>::iterator itB = ids_of_species[species_B].begin(); itB != ids_of_species[species_B].end(); itB++){
                  int jatom = (*itB);
                  float jrx = atoms[jatom].rx;
                  float jry = atoms[jatom].ry;
                  float jrz = atoms[jatom].rz;
                  float dx = jrx - irx;
                  dx -= lx * round(dx * ilx);
                  float dy = jry - iry;
                  dy -= ly * round(dy * ily);
                  float dz = jrz - irz;
                  dz -= lz * round(dz * ilz);
                  float rr2 = dx * dx + dy * dy + dz * dz;
                  if (rr2 < r_max2){
                     if (jatom == iatom) continue;
                     if (interact_flag == 1 && atoms[iatom].mol == atoms[jatom].mol) continue;
                     if (interact_flag == 2 && atoms[iatom].mol != atoms[jatom].mol) continue;
                     int ibin = (int)(sqrt(rr2)*ir_bin);
                     partial_hist_ij[ species_A ][ species_B ][ ibin ] += 1;
                  }
               }
            }
         }
      }
      dump_file.close();
   }

   volume /= (float)n_samples;
   float rho = (float)n_atoms / volume;
   cout<<"\nNumber of samples examined: "<<n_samples<<endl;
   cout<<"Ensemble averaged System volume: "<<volume<<endl;
   cout<<endl;

   // Output the gij's of the selected species
   for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();
        it_pair != pairs.end(); it_pair++){
      int species_A = it_pair->first;
      int species_B = it_pair->second;
      for (int ibin = 0; ibin < n_bin; ibin++){
         float shell_vol = 4.0 / 3.0 * M_PI * (pow(ibin+1, 3) - pow(ibin, 3)) * pow(r_bin, 3);
         float nid = shell_vol * rho;
         float NA = (float)ids_of_species[species_A].size();
         float NB = (float)ids_of_species[species_B].size();
         //partial_gij[species_A][species_B][ibin] = (float)partial_hist_ij[species_A][species_B][ibin] / ((float)n_samples) * (float)n_atoms / ( NA * N_inter[species_A][species_B] * nid);
         partial_gij[species_A][species_B][ibin] = (float)partial_hist_ij[species_A][species_B][ibin] / ((float)n_samples) * (float)n_atoms / (NA * NB * nid);
         // fixit: since there are no selfinteractions in bulk, shouldn't we scale as follows?
         //if (species_A == species_B) partial_gij[species_A][species_B][ibin] *= (float)NB / ((float)NB - 1.0);
      }
   }

   if (export_verbose) {
      for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();
           it_pair != pairs.end(); it_pair++){
         int species_A = it_pair->first;
         int species_B = it_pair->second;
         string gij_file_name;
         gij_file_name = rdffile + "."+species_of_num[species_A]+"-"+species_of_num[species_B]+".dat";
         FILE * pFile;
         pFile = fopen (gij_file_name.c_str(),"w");
         cout<<"writting gij of pair "<<species_A<<" "<<species_B<<" ( "<<species_of_num[species_A]<<" "<<species_of_num[species_B]<<" )"<<endl;
         for (int ibin = 0; ibin < n_bin; ibin++){
            fprintf (pFile, "%10d %10.6f %10.6f\n",ibin, r_bin*((float)ibin+0.5), partial_gij[species_A][species_B][ibin]);
         }
         fclose(pFile);
      }
   }

   // Output the gij's of the species in a single file
   FILE * pFile;
   string file_all_gijs = rdffile + ".dat";
   cout<<"\nWritting all gij's in file: "<<file_all_gijs<<endl;

   pFile = fopen (file_all_gijs.c_str(),"w");

   fprintf (pFile, "Nall = %-20d\n", n_atoms);
   fprintf (pFile, "rho = %-.5f\n", rho);
   fprintf (pFile, "%-20s", "type");
   for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();it_pair != pairs.end(); it_pair++){
      std::string pair;
      pair += species_of_num[it_pair->first]+"_"+species_of_num[it_pair->second];
      fprintf (pFile, "%-20s", pair.c_str());
   }
   fprintf (pFile, "%-20s ", "\nca");
   for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();it_pair != pairs.end(); it_pair++)
      fprintf (pFile, "%-20.5f", c_of_species[it_pair->first]);
   fprintf (pFile, "%-20s ", "\ncb");
   for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();it_pair != pairs.end(); it_pair++)
      fprintf (pFile, "%-20.5f", c_of_species[it_pair->second]);
   fprintf (pFile, "\n");
   for (int ibin = 0; ibin < n_bin; ibin++){
      fprintf (pFile, "%-20.5f", r_bin*((float)ibin+0.5));
      for (std::vector< std::pair <int, int> >::iterator it_pair = pairs.begin();it_pair != pairs.end(); it_pair++)
         fprintf (pFile, "%-20.5f", partial_gij[it_pair->first][it_pair->second][ibin]);
      fprintf (pFile, "\n");
   }
   fclose(pFile);

   cout<<"Done!"<<endl;

   return 0;
}
