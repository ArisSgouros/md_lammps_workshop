 /*
  * MIT License
  *
  * Copyright (c) 2023 ArisSgouros
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
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <vector>

using namespace std;

std::vector<std::string> tokenize(std::string input_string) {
   std::string buf;
   std::stringstream ss(input_string);
   std::vector<std::string> tokens;
   while (ss >> buf)
      tokens.push_back(buf);
   return tokens;
}

int main(int argc, char** argv){

   // Parse the command line arguments
   if (argc != 4){
      cout<<endl;;
      cout<<"The required formats are the following:"<<endl;
      cout<<endl;
      cout<<"./msd \"parameter_file\" \"dump_file\" \"time_step\""<<endl;
      cout<<endl;
      cout<<"exiting.."<<endl;
      return 0;
   }

   std::string param_filename = argv[1];
   std::string dump_filename  = argv[2];
   double      time_step      = stof(argv[3]);

   std::vector <std::string> lines_of_file;
   std::vector <std::string> tokens;
   std::string current_line;

   // Read the parameter file
   cout<<"\nReading the parameters from file: "<<param_filename<<endl;
   ifstream param_file(param_filename.c_str(), ifstream::in);

   int iline;
   int DUMP_COL_ID=0, DUMP_COL_MOL=0, DUMP_COL_RX=0, DUMP_COL_RY=0, DUMP_COL_RZ=0;
   int mols_line_start;

   lines_of_file.clear();
   iline = 0;
   int n_mols=0;
   while (param_file.good()) {

      // A temporary string for the current line of the file.
      getline(param_file, current_line);

      // Add current line to file's array of lines.
      lines_of_file.push_back(current_line);

      if (current_line.find("number of molecules") != std::string::npos){
         n_mols = atoi(current_line.c_str());
         cout<<"number of molecules: "<<n_mols<<endl;
      }
      if (current_line.find("Molecules") != std::string::npos)
         mols_line_start = iline+2;
      iline += 1;
   }
   param_file.close();

   // Assign local IDs to molecules
   std::vector< int > mol_ids;
   std::vector< int > at_mol;
   std::map<int, int> locid_of_mol;
   int snum = 0;
   for (int iline = mols_line_start; iline < mols_line_start + n_mols; iline++){
      int mol_Id = atoi(lines_of_file[iline].c_str());
      mol_ids.push_back( mol_Id );
      locid_of_mol[mol_Id] = snum;
      at_mol.push_back(0);
      snum ++;
   }
   cout<<"\nlist of mols to be examined"<<endl;
   for (std::vector< int >::iterator it_mol = mol_ids.begin();
        it_mol != mol_ids.end(); it_mol++){
      cout<<(*it_mol)<<" "<<locid_of_mol[(*it_mol)]<<endl;
   }

   // read the format of the lammps dump file and the number of atoms per mol
   {
      cout<<"\nGet the format from the lammps dump file: "<<dump_filename<<endl;
      ifstream dump_file(dump_filename.c_str(), ifstream::in);
      getline(dump_file, current_line); // ITEM: TIMESTEP
      getline(dump_file, current_line);
      getline(dump_file, current_line); // ITEM: NUMBER OF ATOMS
      getline(dump_file, current_line);
      int n_atom = atoi(current_line.c_str());
      getline(dump_file, current_line); // ITEM: BOX BOUNDS pp pp pp
      getline(dump_file, current_line);
      getline(dump_file, current_line);
      getline(dump_file, current_line);
      getline(dump_file, current_line); // ITEM: ATOMS ...
      tokens = tokenize(current_line);
      int pos = -2;
      for (std::vector< std::string >::iterator it_token = tokens.begin();
           it_token != tokens.end(); it_token++){
         if ( (*it_token) == "id" ){
            DUMP_COL_ID = pos;
            cout<<"id col: "<<DUMP_COL_ID<<endl;
         }
         if ( (*it_token) == "mol" ){
            DUMP_COL_MOL = pos;
            cout<<"mol col: "<<DUMP_COL_MOL<<endl;
         }
         if ( (*it_token) == "xu" || (*it_token) == "xs" ){
            DUMP_COL_RX = pos;
            cout<<"x col: "<<DUMP_COL_RX<<endl;
         }
         if ( (*it_token) == "yu" || (*it_token) == "ys"  ){
            DUMP_COL_RY = pos;
            cout<<"y col: "<<DUMP_COL_RY<<endl;
         }
         if ( (*it_token) == "zu" || (*it_token) == "zs"  ){
            DUMP_COL_RZ = pos;
            cout<<"z col: "<<DUMP_COL_RZ<<endl;
         }
         pos++;
      }

      for (int ii = 0; ii < n_atom; ii++){
         getline(dump_file, current_line);
         tokens = tokenize(current_line);
         int imol   = stoi(tokens[DUMP_COL_MOL]);
         if (std::find(mol_ids.begin(), mol_ids.end(), imol) != mol_ids.end()){
            at_mol[locid_of_mol[imol]] ++;
         }
      }
   }

   lines_of_file.clear();
   lines_of_file.push_back(current_line);

   // Initialize a vector of dimenions [frame, mol, coord*3]
   vector< vector < vector<double>>> frameMolRcm;

   // read the center-of-mass coordinates from the lammps file
   cout<<"\nread the center-of-mass coordinates from the lammps file: "<<dump_filename<<endl;
   ifstream dump_file(dump_filename.c_str(), ifstream::in);

   int nFrames = 0;
   while (dump_file.good()) {
      getline(dump_file, current_line); // ITEM: TIMESTEP
      if (current_line.empty()) break;
      getline(dump_file, current_line);
      getline(dump_file, current_line); // ITEM: NUMBER OF ATOMS
      getline(dump_file, current_line);
      int n_atom = atoi(current_line.c_str());
      getline(dump_file, current_line); // ITEM: BOX BOUNDS pp pp pp
      getline(dump_file, current_line);
      getline(dump_file, current_line);
      getline(dump_file, current_line);
      getline(dump_file, current_line); // ITEM: ATOMS ...

      vector < vector<double>> molRcm;
      vector<double> pos{0.0, 0.0, 0.0};
      for (int ii = 0; ii < n_mols; ii++)
         molRcm.push_back(pos);


      for (int ii = 0; ii < n_atom; ii++){
         getline(dump_file, current_line);
         tokens = tokenize(current_line);
         int imol   = stoi(tokens[DUMP_COL_MOL]);

	 if (std::find(mol_ids.begin(), mol_ids.end(), imol) != mol_ids.end()){
            float rx   = stof(tokens[DUMP_COL_RX]);
            float ry   = stof(tokens[DUMP_COL_RY]);
            float rz   = stof(tokens[DUMP_COL_RZ]);
            int loc_id = locid_of_mol[imol];
            molRcm[loc_id][0] += rx;
            molRcm[loc_id][1] += ry;
            molRcm[loc_id][2] += rz;
         }
      }
      for (int ii = 0; ii < n_mols; ii++){
         molRcm[ii][0] /= (double)at_mol[ii];
         molRcm[ii][1] /= (double)at_mol[ii];
         molRcm[ii][2] /= (double)at_mol[ii];
      }
      frameMolRcm.push_back(molRcm);
      nFrames++;
   }

   // compute the mean square displacement
   vector< int > msd_count;
   vector< vector< double>> msd_dim;

   vector< double> pos = {0.0, 0.0, 0.0};
   for (int tt=0; tt<nFrames; tt++){
      msd_dim.push_back(pos);
      msd_count.push_back(0);
   }

   for (int imol=0; imol<n_mols; imol++){
      for (int t0=0; t0<nFrames; t0++){
         for (int tt=0; tt<(nFrames-t0); tt++){
            msd_count[t0] += 1;
            vector< double> r0 = { frameMolRcm[tt+t0  ][imol][0] ,frameMolRcm[tt+t0  ][imol][1], frameMolRcm[tt+t0  ][imol][2]};
            vector< double> r1 = { frameMolRcm[tt     ][imol][0] ,frameMolRcm[tt     ][imol][1], frameMolRcm[tt     ][imol][2]};
            msd_dim[t0][0] += pow(r1[0]-r0[0],2);
            msd_dim[t0][1] += pow(r1[1]-r0[1],2);
            msd_dim[t0][2] += pow(r1[2]-r0[2],2);
         }
      }
   }
   for (int tt=0; tt<nFrames; tt++){
      msd_dim[tt][0] /= (double)msd_count[tt];
      msd_dim[tt][1] /= (double)msd_count[tt];
      msd_dim[tt][2] /= (double)msd_count[tt];
   }

   FILE * pFile;
   string file_msd = "o.msd.dat";
   cout<<"\nWritting msd in file: "<<file_msd<<endl;

   pFile = fopen (file_msd.c_str(),"w");

   fprintf (pFile, "%-20s %-20s %-20s %-20s %-20s\n", "time", "x", "y", "z", "total");
   for (int tt=0; tt<nFrames; tt++)
      fprintf (pFile, "%-20.6f %-20.6f %-20.6f %-20.6f %-20.6f\n",time_step*(double)tt, msd_dim[tt][0], msd_dim[tt][1], msd_dim[tt][2], msd_dim[tt][0]+msd_dim[tt][1]+msd_dim[tt][2]);
   fclose(pFile);

   return 0;
}
