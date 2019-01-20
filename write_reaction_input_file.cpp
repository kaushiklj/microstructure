#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void write_reaction_input_file(int nframes1, int natoms1,int number_molecules, atom **at_list1, molecule **molecule_list1){
     int i,j,k;
     int n1,n2,wt;
     int frames;
     string s1,s2;
     ofstream file1;
     
     file1.open("reaction.in",ios::app);
     cout << "Writing reaction input file\n";
     frames=0;
     for(i=0;i<natoms1;i++){
                            n1 = at_list1[frames][i].return_mol_no();
                            n1 = n1 - 1;
                            s1 = molecule_list1[frames][n1].return_molecule_name();
                            wt = molecule_list1[frames][n1].return_molecule_weight();           
                            file1 << setw(2) << at_list1[frames][i].return_atomname() << " " << setw(8) << i << " " 
                                  << setw(2) << at_list1[frames][i].num_neighbors << " " << setw(6) << at_list1[frames][i].return_mol_no() << "  " << setw(50) << s1 << " " << wt << endl;                             
                            }
     file1.close();
     }
